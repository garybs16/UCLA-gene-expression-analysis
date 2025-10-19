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
