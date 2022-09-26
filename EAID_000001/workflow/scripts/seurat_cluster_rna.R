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

proj <- readRDS(file = input_paths[["project_in"]])
print(proj) ####

# ref <- readRDS(file = input_paths[["project_ref"]])
# print(head(ref@meta.data)) ####
# print(ref) ####
# ref <- subset(x = ref, subset = `Factor.Value.inferred.cell.type...authors.labels.` != "")
# print(ref) ####
# print(proj) ####

# anchors <- FindTransferAnchors(
#   reference = ref,
#   query = proj,
#   normalization.method = "SCT",
#   dims = 1:30, 
#   reference.reduction = "pca",
#   reference.assay = "SCT",
#   query.assay = "integrated",
#   reduction = "pcaproject"
# )
# proj <- MapQuery(
#   anchorset = anchors,
#   query = proj,
#   reference = ref,
#   refdata = "cell_type",
#   reference.reduction = "pca"
#   # reduction.model = "umap"
# )
# print(head(proj@meta.data)) ####
# print(head(ref@meta.data)) ####
# proj$cell_type <- proj$predicted.id

# plt_ref <- DimPlot(proj, reduction = "umap.ref", group.by = "cell_type")
# ggsave(output_paths[["umap_ref"]], plt, device = "pdf")

# all.genes <- rownames(proj)
# proj <- ScaleData(proj, features = all.genes)

# proj <- RunPCA(proj, features = VariableFeatures(object = proj))

proj <- FindNeighbors(proj, dims = 1:30, reduction="harmony")
proj <- FindClusters(object = proj) 
# proj <- RunUMAP(proj, dims = 1:30)

plt <- DimPlot(proj, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
ggsave(output_paths[["umap"]], plt, device = "pdf")

# plt <- DimPlot(proj, reduction = "umap", group.by = "dataset")
# ggsave(output_paths[["umap_dataset"]], plt, device = "pdf")

saveRDS(proj, file = output_paths[["project_out"]])