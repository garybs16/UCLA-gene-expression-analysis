# Tahoe-100M Analysis Summary

## Run Snapshot

- Cells analyzed: 20,000
- Embedding dimensions: 50
- Clusters discovered: 10
- Cumulative explained variance across retained components: 0.095

## Cluster Quality Metrics

- Silhouette score: 0.074
- Davies-Bouldin score: 2.579
- Calinski-Harabasz score: 1146.3
- Cluster balance ratio: 0.168

## Cluster Sizes

- Cluster 0: 3,478 cells
- Cluster 1: 1,168 cells
- Cluster 2: 618 cells
- Cluster 3: 851 cells
- Cluster 4: 2,165 cells
- Cluster 5: 1,012 cells
- Cluster 6: 1,849 cells
- Cluster 7: 3,675 cells
- Cluster 8: 3,398 cells
- Cluster 9: 1,786 cells

## Metadata Coverage

- Unique `drug` values: 44
- Unique `cell_line_id` values: 50
- Unique `moa-fine` values: 8
- Unique `plate` values: 1

## Most Common Drugs in Sample

- Quinestrol: 736 cells
- Aliskiren: 650 cells
- Retinoic acid: 637 cells
- Riluzole hydrochloride: 623 cells
- (R)-Verapamil (hydrochloride): 603 cells
- Pasireotide (acetate): 602 cells
- Isocorydine: 564 cells
- Palmatine (chloride): 551 cells
- Berbamine: 547 cells
- Tolcapone: 543 cells

## Interview Talking Points

- The pipeline is memory-aware: it streams cells and keeps the expression matrix sparse.
- The normalization mirrors a common single-cell RNA-seq workflow: library-size scaling plus log1p.
- TruncatedSVD is used because the normalized expression matrix remains sparse.
- The clustering stage creates a compact representation for downstream drug-response exploration.
- Cluster quality metrics are exported so the result can be evaluated beyond visual inspection.
- Current limitations: clustering is unsupervised, results should be validated against biological labels, and a larger run would need batch-effect controls.
