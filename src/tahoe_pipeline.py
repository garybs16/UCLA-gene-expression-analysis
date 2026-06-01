"""Build low-dimensional single-cell embeddings from Tahoe-100M.

The pipeline streams expression profiles from Hugging Face, constructs a sparse
cell-by-gene matrix, applies CPM-style log normalization, reduces dimensions
with TruncatedSVD, and clusters cells in embedding space.
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from datasets import load_dataset
from scipy.sparse import csr_matrix
from sklearn.cluster import KMeans
from sklearn.decomposition import TruncatedSVD


DEFAULT_DATASET = "tahoebio/Tahoe-100M"


@dataclass(frozen=True)
class PipelineConfig:
    dataset: str
    samples: int
    pca_dim: int
    clusters: int
    output_dir: str
    random_state: int


def load_gene_vocab(dataset_name: str) -> dict[int, str]:
    """Load token-to-Ensembl gene metadata."""
    gene_meta = load_dataset(dataset_name, name="gene_metadata", split="train").to_pandas()
    required_cols = {"token_id", "ensembl_id"}
    missing = required_cols.difference(gene_meta.columns)
    if missing:
        raise ValueError(f"gene metadata is missing columns: {sorted(missing)}")

    return {
        int(row["token_id"]): str(row["ensembl_id"])
        for _, row in gene_meta.iterrows()
    }


def build_sparse_expression_matrix(
    cell_stream: Iterable[dict],
    gene_vocab: dict[int, str],
    n_cells: int,
) -> tuple[csr_matrix, pd.DataFrame, pd.DataFrame]:
    """Convert streamed tokenized cells into a sparse cell-by-gene matrix."""
    token_ids, gene_ids = zip(*sorted(gene_vocab.items()))
    token_to_col = {token_id: col for col, token_id in enumerate(token_ids)}

    data: list[float] = []
    indices: list[int] = []
    indptr: list[int] = [0]
    observations: list[dict] = []

    for i, cell in enumerate(cell_stream):
        if i >= n_cells:
            break

        genes = cell.get("genes", []) or []
        expressions = cell.get("expressions", []) or []

        # Some records include a leading sentinel value; it is not expression.
        if expressions and expressions[0] < 0:
            genes = genes[1:]
            expressions = expressions[1:]

        for gene_token, expression in zip(genes, expressions):
            col = token_to_col.get(gene_token)
            if col is not None:
                indices.append(col)
                data.append(float(expression))

        indptr.append(len(data))
        observations.append(
            {key: value for key, value in cell.items() if key not in {"genes", "expressions"}}
        )

    matrix = csr_matrix((data, indices, indptr), shape=(len(indptr) - 1, len(gene_ids)))
    obs_df = pd.DataFrame(observations)
    var_df = pd.DataFrame({"ensembl_id": gene_ids})
    return matrix, obs_df, var_df


def normalize_cpm_log1p(matrix: csr_matrix) -> csr_matrix:
    """Apply counts-per-million normalization followed by log1p."""
    cell_sums = np.asarray(matrix.sum(axis=1)).reshape(-1, 1)
    cell_sums[cell_sums == 0] = 1.0

    normalized = matrix.multiply(1_000_000.0).multiply(1.0 / cell_sums)
    normalized.data = np.log1p(normalized.data)
    return normalized


def run_pipeline(config: PipelineConfig) -> dict:
    """Run the full embedding and clustering workflow."""
    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[info] streaming {config.samples:,} cells from {config.dataset}")
    cell_stream = load_dataset(config.dataset, split="train", streaming=True)

    print("[info] loading gene metadata")
    gene_vocab = load_gene_vocab(config.dataset)

    expression_matrix, obs_df, var_df = build_sparse_expression_matrix(
        cell_stream=cell_stream,
        gene_vocab=gene_vocab,
        n_cells=config.samples,
    )
    print(
        "[info] built sparse matrix: "
        f"cells={expression_matrix.shape[0]:,}, "
        f"genes={expression_matrix.shape[1]:,}, "
        f"nonzeros={expression_matrix.nnz:,}"
    )

    normalized = normalize_cpm_log1p(expression_matrix)

    n_components = max(2, min(config.pca_dim, normalized.shape[1] - 1))
    print(f"[info] reducing to {n_components} dimensions with TruncatedSVD")
    svd = TruncatedSVD(n_components=n_components, random_state=config.random_state)
    embedding = svd.fit_transform(normalized)

    n_clusters = min(config.clusters, max(2, embedding.shape[0] // 50))
    print(f"[info] clustering embeddings with KMeans(k={n_clusters})")
    kmeans = KMeans(
        n_clusters=n_clusters,
        n_init="auto",
        random_state=config.random_state,
    )
    labels = kmeans.fit_predict(embedding)

    embedding_df = pd.DataFrame(
        embedding,
        columns=[f"pc{i + 1}" for i in range(embedding.shape[1])],
    )
    output_df = pd.concat([obs_df.reset_index(drop=True), embedding_df], axis=1)
    output_df["cluster"] = labels

    embeddings_path = output_dir / "embeddings.csv"
    variance_path = output_dir / "pca_variance.csv"
    metadata_path = output_dir / "run_metadata.json"

    output_df.to_csv(embeddings_path, index=False)
    pd.DataFrame(
        {
            "pc": [f"pc{i + 1}" for i in range(len(svd.explained_variance_ratio_))],
            "explained_variance_ratio": svd.explained_variance_ratio_,
        }
    ).to_csv(variance_path, index=False)

    metadata = {
        "config": asdict(config),
        "cells": int(expression_matrix.shape[0]),
        "genes": int(expression_matrix.shape[1]),
        "nonzero_expression_values": int(expression_matrix.nnz),
        "density": float(expression_matrix.nnz / (expression_matrix.shape[0] * expression_matrix.shape[1])),
        "components": int(n_components),
        "clusters": int(n_clusters),
        "cumulative_explained_variance": float(np.sum(svd.explained_variance_ratio_)),
    }
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    print(f"[saved] {embeddings_path}")
    print(f"[saved] {variance_path}")
    print(f"[saved] {metadata_path}")
    return metadata


def parse_args() -> PipelineConfig:
    parser = argparse.ArgumentParser(
        description="Stream Tahoe-100M cells and build clustered transcriptomic embeddings.",
    )
    parser.add_argument("--dataset", default=DEFAULT_DATASET)
    parser.add_argument("--samples", type=int, default=20_000)
    parser.add_argument("--pca-dim", type=int, default=50)
    parser.add_argument("--clusters", type=int, default=10)
    parser.add_argument("--output-dir", default="outputs")
    parser.add_argument("--random-state", type=int, default=42)
    args = parser.parse_args()

    if args.samples < 2:
        raise ValueError("--samples must be at least 2")
    if args.pca_dim < 2:
        raise ValueError("--pca-dim must be at least 2")
    if args.clusters < 2:
        raise ValueError("--clusters must be at least 2")

    return PipelineConfig(
        dataset=args.dataset,
        samples=args.samples,
        pca_dim=args.pca_dim,
        clusters=args.clusters,
        output_dir=args.output_dir,
        random_state=args.random_state,
    )


def main() -> None:
    os.environ.setdefault("HF_HUB_ENABLE_HF_TRANSFER", "1")
    os.environ.setdefault("HF_HUB_READ_TIMEOUT", "60")
    run_pipeline(parse_args())


if __name__ == "__main__":
    main()
