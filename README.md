# 🧬 UCLA Tahoe-100M Gene Expression Analysis System

**A scalable pipeline for streaming, transforming, and clustering massive single-cell gene expression datasets from [Tahoe-100M](https://huggingface.co/datasets/tahoebio/Tahoe-100M).**

---

## 📘 Overview

This repository provides a **lightweight, high-performance analysis system** for the **UCLA Tahoe-100M** dataset — a 100-million-cell benchmark dataset for large-scale gene expression research.

The script `ucla_tahoe100m_pipeline_v2.py` supports:
- Streaming cell data directly from Hugging Face without full dataset download  
- Sparse matrix construction and normalization  
- Dimensionality reduction using **TruncatedSVD (sparse PCA)**  
- **K-Means clustering** in reduced space  
- Automatic export of embeddings and explained variance

This pipeline is ideal for exploratory gene expression analysis, large-scale benchmarking, and downstream bioinformatics workflows.

---

## 🚀 Quick Start

### **1. Install dependencies**
```bash
pip install datasets numpy pandas scipy scikit-learn



2. Run the pipeline
python ucla_tahoe100m_pipeline_v2.py --samples 20000 --pca_dim 50 --clusters 10

3. Arguments
Argument	Default	Description
--dataset	tahoebio/Tahoe-100M	Hugging Face dataset ID
--samples	20000	Number of cells to stream
--pca_dim	50	Number of principal components (SVD dims)
--clusters	10	K-Means cluster count
📊 Outputs

The pipeline produces two CSV files after completion:

File	Description
embeddings.csv	Principal components + metadata + cluster assignment for each cell
pca_variance.csv	Explained variance ratio for each principal component
🧠 Methodology

Streaming and Sparse Matrix Construction
Loads the Tahoe-100M dataset in streaming mode using Hugging Face’s datasets API.
Builds a sparse CSR matrix of expression counts indexed by Ensembl gene IDs.

Normalization
Performs CPM-like normalization per cell and applies a log1p transformation to control scale.

Dimensionality Reduction
Uses TruncatedSVD (a sparse PCA equivalent) to project gene expression data into lower-dimensional space.

Clustering
Applies K-Means clustering to the reduced embedding space to identify putative cell groups or subtypes.

🧩 Example Output

Example of pca_variance.csv:

pc	explained_variance_ratio
pc1	0.123
pc2	0.089
pc3	0.065
...	...
🧪 Dataset Reference

Tahoe-100M — tahoebio/Tahoe-100M on Hugging Face

A 100M-cell dataset for standardized, large-scale gene expression analysis.
Includes tokenized gene vocabularies, metadata, and expression matrices optimized for streaming and scalable computation.

⚙️ System Requirements
Resource	Recommended
CPU	8+ cores
RAM	≥16 GB (scales with --samples)
Disk	~2 GB temporary storage
OS	Linux / macOS / WSL2
📈 Example Workflow
# Stream 50k cells and compute 100 principal components with 20 clusters
python ucla_tahoe100m_pipeline_v2.py --samples 50000 --pca_dim 100 --clusters 20

# Inspect PCA variance
cat pca_variance.csv | head

# Check embeddings and cluster labels
head embeddings.csv

🧰 Citation

If you use this code or dataset in your research:

@dataset{tahoebio_Tahoe_100M,
  title     = {Tahoe-100M: A Large-Scale Gene Expression Benchmark Dataset},
  author    = {TahoeBio and UCLA Computational Biology},
  year      = {2024},
  url       = {https://huggingface.co/datasets/tahoebio/Tahoe-100M}
}
