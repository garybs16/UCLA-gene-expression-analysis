# Interview Notes

## 30-Second Project Pitch

I built a Python pipeline for analyzing single-cell transcriptomic response data from Tahoe-100M. The pipeline streams cells from Hugging Face, constructs a sparse cell-by-gene matrix, applies CPM-style log normalization, reduces dimensionality with TruncatedSVD, and clusters cells with KMeans. I also added reporting code that generates cluster plots, variance diagnostics, and a drug-by-cluster heatmap for interpretation.

I added a supervised benchmark as a sanity check: using only the learned embeddings, the project trains a logistic regression classifier to test whether mechanism-of-action labels can be predicted better than a majority-class baseline.

## Why These Technical Choices

- Streaming avoids downloading the full Tahoe-100M dataset before analysis.
- Sparse matrices are necessary because most genes are not observed in each cell.
- CPM plus log1p is a standard normalization pattern for single-cell count data.
- TruncatedSVD works directly with sparse matrices, unlike standard PCA workflows that often densify data.
- KMeans gives a simple baseline clustering method that is easy to explain and evaluate.
- Cluster quality metrics give a quantitative check beyond visual inspection.

## Current Run

- 20,000 cells analyzed
- 50 embedding dimensions
- 10 discovered clusters
- 44 drugs represented
- 50 cell lines represented
- 8 mechanism-of-action labels represented
- Cluster validation exported with silhouette, Davies-Bouldin, Calinski-Harabasz, and balance metrics
- Supervised mechanism-of-action benchmark included

## Strong Answers to Likely Questions

**What problem does this solve?**

It creates a compact representation of high-dimensional single-cell gene-expression profiles so downstream analysis can compare drug and cell-line response patterns.

**How do you know the clusters are meaningful?**

I treat them as candidate structure, not final truth. The project exports visualizations, cluster-size balance, and standard clustering metrics. The next stronger validation would compare clusters against known mechanism-of-action labels or train a supervised benchmark from the embeddings.

**Why add a supervised benchmark?**

Clustering is useful, but it can be hard to judge alone. The supervised benchmark asks whether the embeddings preserve label-relevant information by predicting mechanism-of-action classes from the low-dimensional representation.

**Why not ordinary PCA?**

The expression matrix is sparse and high-dimensional. TruncatedSVD is a practical PCA-like method that can operate efficiently on sparse matrices.

**How would you improve it with more time?**

I would add biological validation against known cell-line or mechanism-of-action labels, batch-effect correction, configurable feature filtering, and a supervised evaluation task such as predicting mechanism of action from embeddings.

**What is the biggest limitation?**

The clusters are unsupervised, so they should be interpreted as candidate structure rather than biological truth until validated against labels or external knowledge.

**How would this scale?**

The current pipeline already streams input and stores the expression matrix sparsely. For much larger runs, I would process in chunks, persist intermediate sparse matrices, and consider incremental dimensionality reduction or distributed execution.
