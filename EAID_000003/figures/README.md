This folder contains relevant visualizations of the data.

Each sample (W62, W64, W73) has the following quality control figures:
- ATAC_QC: pre-filtered TSS Enrichment versus unique fragments
- filtered_ATAC_QC: post-filtering TSS Enrichment versus unique fragments
- ATAC_frags: pre-filtered fragment size distribution
- RNA_QC: pre-filtered unique RNA counts, RNA counts, percent mitochondrial

The final integrated and labeled dataset has the following figures:
- UMAP: UMAP labeled by cluster
- UMAP_nolegend: UMAP labeled by cluster without legend
- marker_expression: dotplot of marker gene expression by cluster
- cluster_ratios: barplot of cluster proportions across each sample

The integrated and labeled dataset before high ribosomal and putative doublet clusters are removed:
- presubset_UMAP: UMAP labeled by cluster
- presubset_UMAP_nolegend: UMAP labeled by cluster without legend
- presubset_marker_expression: dotplot of marker gene expression by cluster
- presubset_cluster_ratios: barplot of cluster proportions across each sample
