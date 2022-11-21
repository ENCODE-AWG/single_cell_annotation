# ENCODE Manual Analyses: Pancreas RNA & Multiome Datasets

Analyses for integrating ENCODE scRNA and 10x Multiome datasets

Contact information: 
Ben Parks 
bparks@stanford.edu

## Running the analyses

All analysis steps from data download to processed outputs are summarized in the `run_all.sh` script.

Note that this script may take very long to run end-to-end, and assume high memory and CPU count availability.
It is recommended to run the steps one at a time in a supervized manner.

Standardized outputs are saved in the `out` folder.
Figures are saved in the `figures` folder.
Intermediate data objects are stored in the `data` folder.

## Analysis Steps
- Cells are filtered based on fragment count, TSS enrichment, and/or UMI count (parameters in config/sample_cutoffs.tsv)
- Run sctransform followed by Seurat anchor integration for the RNA (code in `analysis/01b_integrate_samples.R`)
- Use a combination of Seurat alignment with outside datasets (Tabula Sapiens, Azimuth) and
  marker gene analysis (via Seurat AddModuleScore) to attach cell types to the unbiased clusters

## Figure descriptions
- `Azimuth_clusters.pdf` Heatmap of jaccard overlap for cell types annotated by Azimuth alinment vs.
   unbiased cluster IDs from initial integration
- `Azimuth_modules.pdf` Dotplot of module scores for Azimuth marker genes vs. unbiased cluter IDs from
   initial integration
- `Hubmap_markers.pdf` Dotplot of gene expression for selected markers from the Hubmap ASCT+B tables vs.
   unbiased cluster IDs from initial integration
- `TS_cluters.pdf` Same as `Azimuth_clusters.pdf`, but with alignment to Tabula Sapiens rather than Azimuth
- `umap_annotated_clusts.png` UMAP plot of the integrated datasets (based on RNA), labeled by annotated cell type
- `umap_automated_clusts.png` UMAP plot of the integrated datasets (based on RNA), labeled by automated cluter ID