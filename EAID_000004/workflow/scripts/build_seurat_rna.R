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

# Load the proj dataset
expression_matrix <- ReadMtx(
  mtx = input_paths[["mat"]], 
  features = input_paths[["features"]],
  cells = input_paths[["cells"]]
)
# Initialize the Seurat object with the raw (non-normalized data).
proj <- CreateSeuratObject(counts = expression_matrix, project = params[["sample_name"]], min.cells = 3, min.features = 0)
proj <- AddMetaData(
  object = proj,
  metadata = rep(params[["sample_name"]], length(Cells(proj))),
  col.name = 'dataset'
)
proj[["percent.mt"]] <- PercentageFeatureSet(proj, pattern = "^MT-")
proj[["percent.ribo"]] <- PercentageFeatureSet(proj, pattern = "^RP[SL]")

plt <- VlnPlot(proj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(output_paths[["qc_violin"]], plt, device = "pdf")

plot1 <- FeatureScatter(proj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(proj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plt <- plot1 + plot2
ggsave(output_paths[["qc_scatter"]], plt, device = "pdf")

write.table(proj@meta.data, file=output_paths[["metadata"]], quote=FALSE, sep='\t', col.names = NA)

proj <- subset(proj, subset = nCount_RNA > params[["min_count_rna"]] & percent.mt < params[["max_pct_mito_rna"]])

proj <- SCTransform(proj, vars.to.regress = "percent.mt", verbose = FALSE)
proj <- RunPCA(proj, features = VariableFeatures(object = proj))
proj <- FindNeighbors(proj, dims = 1:30)
proj <- FindClusters(object = proj) 

saveRDS(proj, file = output_paths[["project_out"]])

sink(type = "message")
sink()