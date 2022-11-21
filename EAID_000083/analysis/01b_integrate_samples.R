suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(SeuratDisk)
  library(patchwork)
  library(data.table)
  library(GenomicRanges)
})
# Prevent BLAS ops from trying to oversubscribe cores
RhpcBLASctl::blas_set_num_threads(2)

meta <- read_tsv("config/panc_samples.tsv")
panc_all <- meta %>% pull(id)

seurat_list <- readRDS("data/seurat_raw.rds")

# About 3 minutes
system.time(
  seurat_list <- mclapply(seurat_list, mc.cores=8, function(proj) {
    proj %>%
      SCTransform(method = "glmGamPoi") %>%
      RunPCA()
  })  
)

# About 1 minute
system.time(
  seurat_list <- mclapply(seurat_list, mc.cores=8, function(proj) {
    proj %>%
      FindNeighbors(dims=1:30) %>%
      RunUMAP(dims=1:30) %>%
      FindClusters()
  })  
)


cell_mapping <- read_tsv("data/tabula_sapiens_cell_labels.tsv") %>%
  mutate(cell_label=as.factor(cell_label))

unique_clusts <- unique(cell_mapping$cell_label)
ts_colors <- set_names(clust_colors[seq_along(unique_clusts)], unique_clusts)

# 10 seconds
system.time(
  seurat_list <- mclapply(seurat_list, mc.cores=8, function(proj) {
    proj[["ts_cell_mapping"]] <- cell_mapping %>% filter(sample == proj$sample[1]) %>%
      pull(cell_label, name=cell_id) %>%
      .[colnames(proj)]
    proj
  })
)

elbow_plots <- map(seq_along(seurat_list), function(i) ElbowPlot(seurat_list[[i]], ndims=50) + ggtitle(panc_all[i])) %>%
  patchwork::wrap_plots()

clust_colors <- c(ArchR:::ArchRPalettes$stallion, ArchR:::ArchRPalettes$kelly)
names(clust_colors) <- NULL

umap_individual_plots <- map(seq_along(seurat_list), function(i) {
  DimPlot(seurat_list[[i]], label=TRUE) + 
    ggtitle(panc_all[i]) +
    scale_color_manual(values=clust_colors) +
    guides(color="none")
}) %>%
  patchwork::wrap_plots()


umap_cell_mapping <- map(seq_along(seurat_list), function(i) {
  unique_clusts <- unique(cell_mapping$cell_label)
  colors <- set_names(clust_colors[seq_along(unique_clusts)], unique_clusts)
  DimPlot(seurat_list[[i]], group.by="ts_cell_mapping") + 
    ggtitle(panc_all[i]) +
    scale_color_manual(values=colors)
}) %>%
  patchwork::wrap_plots() + patchwork::plot_layout(guides="collect")

saveRDS(seurat_list, "data/seurat_sct_list.rds", compress=FALSE)

pdf("plots/sct_per_sample_plots.pdf", width=14, height=8)
elbow_plots
umap_individual_plots
umap_cell_mapping
dev.off()

features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)

#30 minutes
system.time(
  integration.anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", anchor.features = features)
)
# This also takes several minutes, ot sure of exact timing
seurat_combined <- IntegrateData(anchorset = integration.anchors, normalization.method = "SCT")

# 7 min compressed, 25s uncompressed, definitely not worth waiting 6 minutes to save 4.5GB disk space
system.time(
  saveRDS(seurat_combined, "data/seurat_combined.rds", compress=FALSE)
)


# Liberate some threads for PCA
RhpcBLASctl::blas_set_num_threads(8)
# 3 minutes
system.time(
  seurat_combined <- seurat_combined %>%
    RunPCA() %>%
    FindNeighbors(dims=1:30) %>%
    RunUMAP(dims=1:30) %>%
    FindClusters()
)


pdf("plots/combined_summary_plots.pdf", width=14, height=8)

# Some overall UMAP plots & cell type composition
DimPlot(seurat_combined, label = TRUE) + 
  scale_color_manual(values=clust_colors) + 
  guides(color="none") +
  coord_fixed()

DimPlot(seurat_combined, group.by="ts_cell_mapping", label = TRUE) + 
  scale_color_manual(values=ts_colors) + 
  coord_fixed()

DimPlot(seurat_combined, group.by="sample", shuffle=TRUE) + 
  scale_color_brewer(palette="Dark2") + 
  coord_fixed()

clust_mapping <- tibble(
  sample = seurat_combined$sample,
  clust = factor(str_c("C", seurat_combined$seurat_clusters), str_c("C", levels(seurat_combined$seurat_clusters)))
) %>% group_by(sample, clust) %>%
  summarize(count = n()) %>% 
  ungroup()

#cell_type_fraction_plot
clust_mapping %>%
  group_by(sample) %>%
  mutate(fraction=count/sum(count)) %>%
  ggplot(aes(sample, fraction, fill=clust)) +
  geom_col() + 
  scale_fill_manual(values=clust_colors) +
  theme_classic() + 
  coord_flip()

# sample_id_fraction_plot 
clust_mapping %>%
  group_by(clust) %>%
  mutate(fraction=count/sum(count)) %>%
  ggplot(aes(clust, fraction, fill=sample)) +
  geom_col() + 
  scale_fill_brewer(palette="Set1") +
  theme_classic() +
  coord_flip()

dev.off()
