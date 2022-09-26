Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

proj <- readRDS(file = input_paths[["project_rna"]])

ref <- readRDS(file = input_paths[["project_ref"]])
# print(head(ref@meta.data)) ####
# print(ref) ####
# ref <- subset(x = ref, subset = `Factor.Value.inferred.cell.type...authors.labels.` != "")
# print(ref) ####
# print(proj) ####

anchors <- FindTransferAnchors(
  reference = ref,
  query = proj,
  normalization.method = "SCT",
  dims = 1:30, 
  reference.reduction = "pca",
  reduction = "pcaproject"
)
proj <- MapQuery(
  anchorset = anchors,
  query = proj,
  reference = ref,
  refdata = "cell_type",
  reference.reduction = "pcaproject"
)
proj$cell_type_ref <- proj$predicted.id

proj <- FindNeighbors(proj, dims = 1:30, reduction = "pca")
proj <- RunUMAP(proj, dims = 1:30, reduction = "pca")

plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type_ref")
ggsave(output_paths[["umap"]], plt, device = "pdf")

saveRDS(proj, file = output_paths[["project_out"]])