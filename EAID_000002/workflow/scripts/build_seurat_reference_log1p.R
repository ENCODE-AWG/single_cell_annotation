Sys.setenv(CONDA_BUILD_SYSROOT = "/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

params <- snakemake@params
input_paths <- snakemake@input
output_paths <- snakemake@output
log_paths <- snakemake@log

set.seed(params[["seed"]])

# Load the proj dataset
expression_matrix <- ReadMtx(
  mtx = input_paths[["mat"]], 
  features = input_paths[["features"]],
  cells = input_paths[["cells"]],
  feature.column = 1,
)
# print(expression_matrix) ####
metadata <- read.table(file = input_paths[["metadata"]], sep = ',', header = TRUE)
rownames(metadata) <- metadata$cell
print(head(metadata)) ####
# expression_matrix <- expression_matrix[rownames(metadata)]

# Initialize the Seurat object with the raw (non-normalized data).
proj <- CreateSeuratObject(
    counts = expression_matrix, 
    project = "reference", 
    min.cells = 3, 
    min.features = 200, 
    meta.data = metadata
)
proj[["percent.mt"]] <- PercentageFeatureSet(proj, pattern = "^MT-")
proj$cell_type <- proj[["annot"]]
# proj$cell_type <- replace(proj$cell_type, proj$cell_type == "", "Unknown")
# proj <- subset(proj, subset = cell_type != "NA")
proj <- subset(proj, subset = typeSample == "nucSeq")

plt <- VlnPlot(proj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(output_paths[["qc_violin"]], plt, device = "pdf")

plot1 <- FeatureScatter(proj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(proj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plt <- plot1 + plot2
ggsave(output_paths[["qc_scatter"]], plt, device = "pdf")

proj <- subset(proj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

proj <- NormalizeData(proj, normalization.method = "LogNormalize", scale.factor = 10000)
proj <- FindVariableFeatures(proj, selection.method = "vst", nfeatures = 2000)
proj <- ScaleData(proj)

proj <- RunPCA(proj, features = VariableFeatures(object = proj))

proj <- FindNeighbors(proj, dims = 1:30)

proj <- RunUMAP(proj, dims = 1:30, return.model = TRUE)

plt <- DimPlot(proj, reduction = "umap", group.by = "cell_type")
ggsave(output_paths[["umap"]], plt, device = "pdf")

saveRDS(proj, file = output_paths[["project_out"]])

sink(type = "message")
sink()