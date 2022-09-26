Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])


projs <- lapply(input_paths[["projects_in"]], readRDS)
# print(projs[[1]]) ####
lapply(projs, print) ####

# features <- SelectIntegrationFeatures(object.list = projs, assay = rep("SCT", times = length(projs)))
# projs <- PrepSCTIntegration(object.list = projs, anchor.features = features, assay = "SCT")
# anchors <- FindIntegrationAnchors(object.list = projs, anchor.features = features, reduction = "rpca", normalization.method = "SCT", dims = 1:30, k.anchor = 30)
# proj_merged <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30) 
# DefaultAssay(proj_merged) <- "integrated"

proj_merged <- merge(projs[[1]], projs[-1], project = "merged_rna", add.cell.ids = unlist(params[["samples"]]))

proj_merged <- SCTransform(proj_merged, vars.to.regress = "percent.mt", verbose = FALSE)
proj_merged <- RunPCA(proj_merged)

proj_merged <- FindNeighbors(proj_merged, dims = 1:30, reduction = "pca")
proj_merged <- RunUMAP(proj_merged, dims = 1:30, reduction = "pca")

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_ref")
ggsave(output_paths[["umap_pre_harmony"]], plt, device = "pdf")

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "dataset")
ggsave(output_paths[["umap_dataset_pre_harmony"]], plt, device = "pdf")

# plt <- DimPlot(proj_merged, reduction = "pca", group.by = "dataset")
# ggsave(output_paths[["pca_pre_harmony"]], plt, device = "pdf")

proj_merged <- RunHarmony(proj_merged, "dataset", assay.use = "SCT")

# plt <- DimPlot(proj_merged, reduction = "harmony", group.by = "dataset")
# ggsave(output_paths[["pca_post_harmony"]], plt, device = "pdf")

proj_merged <- FindNeighbors(proj_merged, dims = 1:30, reduction = "harmony")
proj_merged <- RunUMAP(proj_merged, dims = 1:30, reduction = "harmony")

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "cell_type_ref")
ggsave(output_paths[["umap"]], plt, device = "pdf")

plt <- DimPlot(proj_merged, reduction = "umap", group.by = "dataset")
ggsave(output_paths[["umap_dataset"]], plt, device = "pdf")

saveRDS(proj_merged, file = output_paths[["project_out"]])

sink(type = "message")
sink()