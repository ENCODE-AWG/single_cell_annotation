Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(plyr)
library(pheatmap)

confusionMatrix <- function(
  i = NULL, 
  j = NULL
  ){
  ui <- unique(i)
  uj <- unique(j)
  m <- Matrix::sparseMatrix(
    i = match(i, ui),
    j = match(j, uj),
    x = rep(1, length(i)),
    dims = c(length(ui), length(uj))
  )
  rownames(m) <- ui
  colnames(m) <- uj
  m
}

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

proj <- readRDS(file = input_paths[["project_in"]])

cM <- confusionMatrix(proj$seurat_clusters, proj$cell_type_ref)
cM <- as.matrix(cM)
print(cM) ####
print(colnames(cM)) ####
# cM <- cM[, !is.na(colnames(cM))]
cM <- cbind(cM, Unknown = 0.01)
print(cM) ####
labelOld <- rownames(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
proj$labels_named <- plyr::mapvalues(x = proj$seurat_clusters, from = labelOld, to = labelNew)

plt <- pheatmap::pheatmap(
    mat = cM / rowSums(cM), 
    # color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)
ggsave(output_paths[["mat"]], plt)

# new.cluster.ids <- unlist(params[["cluster_names"]])
# names(new.cluster.ids) <- levels(proj)
# proj <- RenameIdents(proj, new.cluster.ids)

plt <- DimPlot(proj, reduction = "umap", label = TRUE, group.by = "labels_named")
ggsave(output_paths[["umap"]], plt)

# proj$labels_named <- Idents(proj)

proj$log10_count_RNA <- log10(proj$nCount_RNA)

plt <- FeaturePlot(object = proj, features = "log10_count_RNA")
ggsave(output_paths[["umap_qc"]], plt)

write.table(proj@meta.data, file=output_paths[["metadata"]], quote=FALSE, sep='\t', col.names = NA)

saveRDS(proj, file = output_paths[["project_out"]])