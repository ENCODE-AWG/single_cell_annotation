suppressPackageStartupMessages({
  library(Seurat)
  library(Azimuth)
  library(SeuratData)
  library(patchwork)
  library(tidyverse)
  library(parallel)
})

annotation <- c(
  "0" = "acinar",
  "1" = "acinar",
  "3" = "acinar",
  "4" = "acinar",
  "5" = "acinar",
  "7" = "acinar",
  "11" = "acinar",
  "12" = "acinar",
  "16" = "acinar",
  "2" = "ductal",
  "6" = "ductal",
  "8" = "ductal",
  "17" = "ductal",
  "18" = "ductal",
  "13" = "beta",
  "22" = "gamma",
  "14" = "alpha_gamma_epsilon",
  "19" = "t",
  "10" = "myeloid",
  "20" = "immune",
  "9" = "stellate",
  "21" = "stellate",
  "25" = "stellate",
  "15" = "endothelial",
  "23" = "endothelial",
  "24" = "schwann"
)
# Load datasets + metadata
seurat_combined <- readRDS("data/seurat_combined.rds")
azimuth_panc <- readRDS("data/azimuth_references/azimuth_panc/ref.Rds")

seurat_combined$annotated_clusters <- annotation[as.character(seurat_combined@active.ident)]
seurat_cells <- str_c(seurat_combined$orig.ident, "#", colnames(seurat_combined)) %>%
  str_remove("_[0-9]$")

annotation_comment <- c(
  "# cell_id: Cell ID used in integrated analysis",
  "# cell_type_id: ENCODE ID of the corresponding cell type object",
  "# cell_type_name: Common name of the cell type",
  "# membership_score: A numeric score for the labeling (Placeholder)"
)
annotations_table <- tibble(
  cell_id = seurat_cells,
  cell_type_id = seurat_combined$annotated_clusters,
  cell_type_name = seurat_combined$annotated_clusters,
  membership_score = NA
)

dir.create("out/labels")
write_lines(annotation_comment, "out/labels/cell_types.tsv.gz")
write_tsv(annotations_table, "out/labels/cell_types.tsv.gz", append=TRUE, col_names=TRUE)

automated_clusters <- tibble(
  cell_id = seurat_cells,
  cell_type_id = seurat_combined@active.ident,
  cell_type_name = seurat_combined@active.ident,
  membership_score = NA
) 
write_lines(annotation_comment, "out/labels/automated_clusters.tsv.gz")
write_tsv(annotations_table, "out/labels/automated_clusters.tsv.gz", append=TRUE, col_names=TRUE)


passing_cells <- fread("out/metadata.tsv.gz") %>%
  as_tibble() %>%
  filter(passed_filtering, !is.na(rna_barcode)) %>%
  pull(cell_id)

meta <- read_tsv("config/panc_samples.tsv")
cutoffs <-read_tsv("config/sample_cutoffs.tsv")

panc_all <- meta$id
panc_multiome <- meta %>% filter(multiome_series != "n/a") %>% pull(id)

gene_metadata <- read_tsv("data/ENCODE/samples/W61_PANC/rna/features.tsv.gz", col_names = c("id", "symbol", "type"))
ensg_to_symbol <- pull(gene_metadata, symbol, name="id")

seurat_mats <- mclapply(panc_all, mc.cores=6, function(s) {
  mat <- readRDS(sprintf("data/seurat/%s/rna_raw.rds", s))
  rownames(mat) <- ensg_to_symbol[rownames(mat)]
  colnames(mat) <- str_c(s, "#", colnames(mat))
  cells <- intersect(colnames(mat), passing_cells)
  mat <- mat[rownames(azimuth_panc),cells]
  mat
})

names(seurat_mats) <- panc_all

# Perform alignment with the azimuth reference
#142 seconds with 4 cores
system.time({
  azimuth_results <- mclapply(seurat_mats, mc.cores=4, function(m) {
    RunAzimuth(CreateSeuratObject(m), "data/azimuth_references/azimuth_panc/")
  })
})

azimuth_cell_type <- lapply(azimuth_results, function(x) x$predicted.annotation.l1) %>%
  set_names(NULL) %>%
  do.call(c, .)


#' Given two lists of alternate cluster assignements (e.g. one integer id per cell), plots a heatmap of the
#' pairwise overlap between clusters as measured by jaccard similarity.
cluster_overlap_heatmap <- function(clust1, clust2) {
  jaccard_similarities <- tidyr::crossing(clust1 = clust1,
                                          clust2 = clust2) %>%
    rowwise() %>%
    mutate(
      jaccard = sum(.env$clust1 == clust1 & .env$clust2 == clust2) / 
        sum(.env$clust1 == clust1 | .env$clust2 == clust2),
      jaccard = replace_na(jaccard, 0)
    )
  
  ggplot(jaccard_similarities, aes(clust1, clust2, fill=jaccard)) +
    geom_tile() +
    scale_fill_gradient(low="white", high="steelblue")
}


# Plot heatmaps for annotation overlaps

dotplot_order <- c("0",  "1",  "3",  "4",  "5",  "7", "11", "12", "16", "2", "6", "8", "17", 
                   "18", "13", "22", "14", "19", "10", "20", "9", "21", "25", 
                     "15", "23",  "24")

azimuth_order <- c("acinar_11", "ductal_10", "beta_12", "alpha_13", "delta_9", "gamma_8", 
                   "epsilon_7", "immune_2", "activated_stellate_6", "quiescent_stellate_3", 
                   "endothelial_4", "schwann_5", "cycling_1")

ts_order <- c("pancreatic acinar cell", "pancreatic ductal cell", "pancreatic beta cell", 
              "pancreatic alpha cell", "pancreatic delta cell", "t cell", "myeloid cell",
              "b cell", "mast cell", "plasma cell", "pancreatic stellate cell", "fibroblast",
              "endothelial cell")

saveRDS(azimuth_cell_type, "data/azimuth_cell_type.rds")

p1 <- cluster_overlap_heatmap(seurat_combined$seurat_clusters, azimuth_cell_type) +
  scale_x_discrete(limits=dotplot_order) +
  scale_y_discrete(limits=str_remove(azimuth_order, "_[0-9]+$")) +
  scale_fill_distiller(palette="Spectral") +
  ggtitle("Azimuth cell calls")

p2 <- cluster_overlap_heatmap(seurat_combined$seurat_clusters, seurat_combined$ts_cell_mapping) +
  scale_x_discrete(limits=dotplot_order) +
  scale_y_discrete(limits=ts_order) +
  scale_fill_distiller(palette="Spectral") +
  ggtitle("TS cell calls")

dir.create("figures")
ggsave(plot = p1, "figures/Azimuth_clusters.pdf", width=7, height=4)
ggsave(plot = p2, "figures/TS_clusters.pdf", width=7, height=4)


# Plot external marker gene lists
azimuth_markers <- read_tsv("config/azimuth_panc_markers.tsv", col_types="")


azimuth_markers_list <- str_split(azimuth_markers$Markers, ", ") %>% 
  set_names(azimuth_markers$Label)

system.time({
  seurat_combined_module <- AddModuleScore(
    seurat_combined, azimuth_markers_list, search = FALSE, name=paste0(names(azimuth_markers_list), "_"))
})


p3 <- DotPlot(seurat_combined_module, cluster.idents=TRUE, features=azimuth_order) +
  coord_flip() +
  theme_classic() +
  scale_y_discrete(limits=dotplot_order) +
  ggtitle("Azimuth marker module scores")
ggsave(plot = p3, "figures/Azimuth_modules.pdf", width=7, height=4)

marker_genes_short <- c(
  "KRT8", "PRSS1", "CFTR", "KRT19", "SYP", "INS", "GCG", "SST", "PTPRC", "CD3E", "CD8A","NKG7", "CD163", "ITGAX",
  "ACTA2", "LAMA2", "PDGFRB", "FABP4", "LEP"
)

p4 <- DotPlot(seurat_combined, cluster.idents=TRUE, 
        features=marker_genes_short) +
  coord_flip() +
  theme_classic() +
  scale_y_discrete(limits=dotplot_order) +
  ggtitle("Hubmap selected marker genes")
ggsave(plot = p4, "figures/Hubmap_markers.pdf", width=7, height=4)


# Plot UMAP

umap_coords <- read_tsv("out/embeddings/umap.tsv.gz", comment="# ")

clust_colors <- c("#D51F26" ,"#272E6A" ,"#208A42" ,"#89288F" ,"#F47D2B" ,"#FEE500" ,"#8A9FD1" ,"#C06CAB" ,"#E6C2DC" ,"#90D5E4" ,"#89C75F"
                  ,"#F37B7D" ,"#9983BD" ,"#D24B27" ,"#3BBCA8" ,"#6E4B9E" ,"#0C727C" ,"#7E1416" ,"#D8A767" ,"#3D3D3D" ,"#FFB300" ,"#803E75"
                  ,"#FF6800" ,"#A6BDD7" ,"#C10020" ,"#CEA262" ,"#817066" ,"#007D34" ,"#F6768E" ,"#00538A" ,"#FF7A5C" ,"#53377A" ,"#FF8E00"
                  ,"#B32851" ,"#F4C800" ,"#7F180D" ,"#93AA00" ,"#593315" ,"#F13A13" ,"#232C16")

p5 <- DimPlot(seurat_combined, label=TRUE) +
  scale_color_manual(values=clust_colors) +
  guides(color="none") +
  coord_fixed()

ggsave(plot=p5, "figures/umap_automated_clusts.png", width = 4, height = 4, dpi=400)

p6 <- DimPlot(seurat_combined, label=TRUE, group.by="annotated_clusters") +
  scale_color_manual(values=clust_colors) +
  guides(color="none") +
  coord_fixed() +
  labs(title = "automated clusters")

ggsave(plot=p6, "figures/umap_annotated_clusts.png", width = 4, height = 4, dpi=400)

