# pip install datasets scipy pandas numpy scikit-learn
import os, argparse
from datasets import load_dataset
import numpy as np, pandas as pd
from scipy.sparse import csr_matrix
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import KMeans

# Make HF downloads more tolerant
os.environ.setdefault("HF_HUB_ENABLE_HF_TRANSFER", "1")
os.environ.setdefault("HF_HUB_READ_TIMEOUT", "60")

def build_matrix(stream, gene_vocab, n_cells):
    # token_id -> col index
    token_ids, gene_ids = zip(*sorted(gene_vocab.items()))
    tok2col = {tid: i for i, tid in enumerate(token_ids)}

    data, indices, indptr, obs = [], [], [0], []
    for i, cell in enumerate(stream):
        if i >= n_cells: break
        genes = cell.get("genes", []) or []
        exprs = cell.get("expressions", []) or []

        # drop possible sentinel/negative first value
        if exprs and exprs[0] < 0:
            genes, exprs = genes[1:], exprs[1:]

        cols = []
        vals = []
        for g, v in zip(genes, exprs):
            if g in tok2col:
                cols.append(tok2col[g])
                vals.append(v)

        data.extend(vals)
        indices.extend(cols)
        indptr.append(len(data))
        obs.append({k: v for k, v in cell.items() if k not in ("genes", "expressions")})

    X = csr_matrix((data, indices, indptr), shape=(len(indptr) - 1, len(gene_ids)))
    obs_df = pd.DataFrame(obs)
    var_df = pd.DataFrame({"ensembl_id": gene_ids})
    return X, obs_df, var_df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="tahoebio/Tahoe-100M")
    ap.add_argument("--samples", type=int, default=20000)
    ap.add_argument("--pca_dim", type=int, default=50)
    ap.add_argument("--clusters", type=int, default=10)
    args = ap.parse_args()

    print(f"[info] streaming {args.samples} cells …")
    ds = load_dataset(args.dataset, split="train", streaming=True)

    print("[info] loading gene metadata …")
    gene_meta = load_dataset(args.dataset, name="gene_metadata", split="train").to_pandas()
    # Expect columns: token_id, ensembl_id
    gene_vocab = {int(row["token_id"]): row["ensembl_id"] for _, row in gene_meta.iterrows()}

    X, obs_df, var_df = build_matrix(ds, gene_vocab, args.samples)
    print(f"[info] built matrix: cells={X.shape[0]} genes={X.shape[1]} (nnz={X.nnz})")

    # CPM-like normalize per cell (row), then log1p; keep sparse
    cell_sums = np.asarray(X.sum(axis=1)).reshape(-1, 1)
    cell_sums[cell_sums == 0] = 1.0
    X_cpm = X.multiply(1_000_000.0).multiply(1.0 / cell_sums)
    X_cpm.data = np.log1p(X_cpm.data)

    # Dimensionality reduction on sparse with TruncatedSVD (PCA for sparse)
    k = max(2, min(args.pca_dim, X_cpm.shape[1] - 1))
    print(f"[info] TruncatedSVD to {k} dims …")
    svd = TruncatedSVD(n_components=k, random_state=42)
    Z = svd.fit_transform(X_cpm)
    var_ratio = svd.explained_variance_ratio_

    # KMeans clustering in embedding space
    k_clusters = min(args.clusters, max(2, Z.shape[0] // 50))
    print(f"[info] KMeans k={k_clusters} …")
    km = KMeans(n_clusters=k_clusters, n_init="auto", random_state=42)
    labels = km.fit_predict(Z)

    # Save outputs
    emb = pd.DataFrame(Z, columns=[f"pc{i+1}" for i in range(Z.shape[1])])
    out = pd.concat([obs_df.reset_index(drop=True), emb], axis=1)
    out["cluster"] = labels
    out.to_csv("embeddings.csv", index=False)
    pd.DataFrame(
        {"pc": [f"pc{i+1}" for i in range(len(var_ratio))],
         "explained_variance_ratio": var_ratio}
    ).to_csv("pca_variance.csv", index=False)

    print("[saved] embeddings.csv, pca_variance.csv")
    print("[done] You can open embeddings.csv to see PCs + clusters per cell.")

if __name__ == "__main__":
    main()
