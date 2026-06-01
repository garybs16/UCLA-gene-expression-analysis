"""Generate figures and a concise analysis report from pipeline outputs."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, silhouette_score


def pc_columns(columns: pd.Index) -> list[str]:
    return sorted(
        [col for col in columns if re.fullmatch(r"pc\d+", col)],
        key=lambda col: int(col[2:]),
    )


def save_explained_variance(variance: pd.DataFrame, figures_dir: Path) -> None:
    x = np.arange(1, len(variance) + 1)
    y = variance["explained_variance_ratio"].to_numpy()
    cumulative = np.cumsum(y)

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(x, y, color="#3B82F6", label="Individual component")
    ax.plot(x, cumulative, color="#111827", marker="o", linewidth=1.5, label="Cumulative")
    ax.set_title("Explained Variance by TruncatedSVD Component")
    ax.set_xlabel("Component")
    ax.set_ylabel("Explained variance ratio")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(figures_dir / "explained_variance.png", dpi=180)
    plt.close(fig)


def save_pc_scatter(embeddings: pd.DataFrame, figures_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    scatter = ax.scatter(
        embeddings["pc1"],
        embeddings["pc2"],
        c=embeddings["cluster"],
        cmap="tab10",
        s=5,
        alpha=0.55,
        linewidths=0,
    )
    ax.set_title("Cell Embeddings Colored by KMeans Cluster")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    legend = ax.legend(*scatter.legend_elements(), title="Cluster", frameon=False, loc="best")
    ax.add_artist(legend)
    ax.grid(alpha=0.18)
    fig.tight_layout()
    fig.savefig(figures_dir / "pc1_pc2_clusters.png", dpi=180)
    plt.close(fig)


def save_cluster_counts(embeddings: pd.DataFrame, figures_dir: Path, tables_dir: Path) -> pd.Series:
    counts = embeddings["cluster"].value_counts().sort_index()
    counts.rename("cell_count").to_csv(tables_dir / "cluster_counts.csv", header=True)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(counts.index.astype(str), counts.values, color="#10B981")
    ax.set_title("Cells per Cluster")
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Cell count")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(figures_dir / "cluster_counts.png", dpi=180)
    plt.close(fig)
    return counts


def save_drug_cluster_heatmap(
    embeddings: pd.DataFrame,
    figures_dir: Path,
    tables_dir: Path,
) -> pd.DataFrame | None:
    if "drug" not in embeddings.columns:
        return None

    top_drugs = embeddings["drug"].value_counts().head(12).index
    subset = embeddings[embeddings["drug"].isin(top_drugs)]
    counts = pd.crosstab(subset["drug"], subset["cluster"])
    proportions = counts.div(counts.sum(axis=1), axis=0)
    counts.to_csv(tables_dir / "drug_cluster_counts.csv")
    proportions.to_csv(tables_dir / "drug_cluster_proportions.csv")

    fig, ax = plt.subplots(figsize=(9, 6))
    image = ax.imshow(proportions.to_numpy(), aspect="auto", cmap="viridis")
    ax.set_title("Cluster Distribution for Most Common Drugs")
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Drug")
    ax.set_xticks(np.arange(len(proportions.columns)), labels=proportions.columns)
    ax.set_yticks(np.arange(len(proportions.index)), labels=proportions.index)
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label("Share of drug's cells")
    fig.tight_layout()
    fig.savefig(figures_dir / "drug_cluster_heatmap.png", dpi=180)
    plt.close(fig)
    return counts


def compute_cluster_metrics(embeddings: pd.DataFrame, tables_dir: Path) -> dict[str, float]:
    pcs = pc_columns(embeddings.columns)
    labels = embeddings["cluster"].to_numpy()
    features = embeddings[pcs].to_numpy()
    n_clusters = len(np.unique(labels))

    if n_clusters < 2 or n_clusters >= len(embeddings):
        metrics = {
            "silhouette_score": np.nan,
            "davies_bouldin_score": np.nan,
            "calinski_harabasz_score": np.nan,
            "cluster_balance_ratio": np.nan,
        }
    else:
        sample_size = min(5_000, len(embeddings))
        metrics = {
            "silhouette_score": float(
                silhouette_score(
                    features,
                    labels,
                    sample_size=sample_size,
                    random_state=42,
                )
            ),
            "davies_bouldin_score": float(davies_bouldin_score(features, labels)),
            "calinski_harabasz_score": float(calinski_harabasz_score(features, labels)),
            "cluster_balance_ratio": float(
                embeddings["cluster"].value_counts().min()
                / embeddings["cluster"].value_counts().max()
            ),
        }

    pd.DataFrame(
        [{"metric": metric, "value": value} for metric, value in metrics.items()]
    ).to_csv(tables_dir / "cluster_quality_metrics.csv", index=False)
    return metrics


def save_metadata_summary(embeddings: pd.DataFrame, tables_dir: Path) -> None:
    fields = [field for field in ["drug", "cell_line_id", "moa-fine", "plate"] if field in embeddings.columns]
    rows = []
    for field in fields:
        rows.append(
            {
                "field": field,
                "unique_values": embeddings[field].nunique(),
                "most_common_value": embeddings[field].value_counts().index[0],
                "most_common_count": int(embeddings[field].value_counts().iloc[0]),
            }
        )
    if rows:
        pd.DataFrame(rows).to_csv(tables_dir / "metadata_summary.csv", index=False)


def write_summary(
    embeddings: pd.DataFrame,
    variance: pd.DataFrame,
    cluster_counts: pd.Series,
    cluster_metrics: dict[str, float],
    report_path: Path,
) -> None:
    pcs = pc_columns(embeddings.columns)
    cumulative_variance = variance["explained_variance_ratio"].sum()
    density_fields = [
        field for field in ["drug", "cell_line_id", "moa-fine", "plate"] if field in embeddings.columns
    ]

    lines = [
        "# Tahoe-100M Analysis Summary",
        "",
        "## Run Snapshot",
        "",
        f"- Cells analyzed: {len(embeddings):,}",
        f"- Embedding dimensions: {len(pcs)}",
        f"- Clusters discovered: {embeddings['cluster'].nunique()}",
        f"- Cumulative explained variance across retained components: {cumulative_variance:.3f}",
        "",
        "## Cluster Quality Metrics",
        "",
        f"- Silhouette score: {cluster_metrics['silhouette_score']:.3f}",
        f"- Davies-Bouldin score: {cluster_metrics['davies_bouldin_score']:.3f}",
        f"- Calinski-Harabasz score: {cluster_metrics['calinski_harabasz_score']:.1f}",
        f"- Cluster balance ratio: {cluster_metrics['cluster_balance_ratio']:.3f}",
        "",
        "## Cluster Sizes",
        "",
    ]
    lines.extend(f"- Cluster {cluster}: {count:,} cells" for cluster, count in cluster_counts.items())

    if density_fields:
        lines.extend(["", "## Metadata Coverage", ""])
        for field in density_fields:
            lines.append(f"- Unique `{field}` values: {embeddings[field].nunique():,}")

    if "drug" in embeddings.columns:
        lines.extend(["", "## Most Common Drugs in Sample", ""])
        for drug, count in embeddings["drug"].value_counts().head(10).items():
            lines.append(f"- {drug}: {count:,} cells")

    lines.extend(
        [
            "",
            "## Interview Talking Points",
            "",
            "- The pipeline is memory-aware: it streams cells and keeps the expression matrix sparse.",
            "- The normalization mirrors a common single-cell RNA-seq workflow: library-size scaling plus log1p.",
            "- TruncatedSVD is used because the normalized expression matrix remains sparse.",
            "- The clustering stage creates a compact representation for downstream drug-response exploration.",
            "- Cluster quality metrics are exported so the result can be evaluated beyond visual inspection.",
            "- Current limitations: clustering is unsupervised, results should be validated against biological labels, and a larger run would need batch-effect controls.",
            "",
        ]
    )

    report_path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create figures and summary from embeddings.csv.")
    parser.add_argument("--embeddings", default="embeddings.csv")
    parser.add_argument("--variance", default="pca_variance.csv")
    parser.add_argument("--output-dir", default="outputs")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    figures_dir = output_dir / "figures"
    tables_dir = output_dir / "tables"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    embeddings = pd.read_csv(args.embeddings)
    variance = pd.read_csv(args.variance)

    required = {"pc1", "pc2", "cluster"}
    missing = required.difference(embeddings.columns)
    if missing:
        raise ValueError(f"embeddings file is missing columns: {sorted(missing)}")

    save_explained_variance(variance, figures_dir)
    save_pc_scatter(embeddings, figures_dir)
    cluster_counts = save_cluster_counts(embeddings, figures_dir, tables_dir)
    save_drug_cluster_heatmap(embeddings, figures_dir, tables_dir)
    cluster_metrics = compute_cluster_metrics(embeddings, tables_dir)
    save_metadata_summary(embeddings, tables_dir)
    write_summary(
        embeddings,
        variance,
        cluster_counts,
        cluster_metrics,
        output_dir / "analysis_summary.md",
    )

    print(f"[saved] {figures_dir / 'explained_variance.png'}")
    print(f"[saved] {figures_dir / 'pc1_pc2_clusters.png'}")
    print(f"[saved] {figures_dir / 'cluster_counts.png'}")
    print(f"[saved] {figures_dir / 'drug_cluster_heatmap.png'}")
    print(f"[saved] {tables_dir / 'cluster_quality_metrics.csv'}")
    print(f"[saved] {tables_dir / 'cluster_counts.csv'}")
    print(f"[saved] {tables_dir / 'metadata_summary.csv'}")
    print(f"[saved] {output_dir / 'analysis_summary.md'}")


if __name__ == "__main__":
    main()
