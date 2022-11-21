suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(dtplyr)
  library(parallel)
  library(Matrix)
  library(Seurat)
})

seurat_combined <- readRDS("data/seurat_combined.rds")
joint_qc <- fread("out/metadata.tsv.gz")

seurat_cells <- str_c(seurat_combined$orig.ident, "#", colnames(seurat_combined)) %>%
  str_remove("_[0-9]$")

all.equal(
  rownames(seurat_combined@reductions$pca@cell.embeddings),
  Cells(seurat_combined))

pca_comment <- c(
  "# cell_id: unique cell identifier for analysis",
  "# PC_*: PCA coordinate, ordered by decreasing variance explained"
)
pca_coords <- as_tibble(seurat_combined@reductions$pca@cell.embeddings) %>%
  mutate(cell_id = seurat_cells) %>%
  select(cell_id, everything())

dir.create("out/embeddings")
write_lines(pca_comment, "out/embeddings/pca.tsv.gz")
write_tsv(pca_coords, "out/embeddings/pca.tsv.gz", append=TRUE, col_names=TRUE)


all.equal(
  rownames(seurat_combined@reductions$umap@cell.embeddings),
  Cells(seurat_combined))
umap_comment <- c(
  "# cell_id: unique cell identifier for analysis",
  "# UMAP_*: UMAP coordinate"
)
umap_coords <- as_tibble(seurat_combined@reductions$umap@cell.embeddings) %>%
  mutate(cell_id = seurat_cells) %>%
  select(cell_id, everything())

write_lines(umap_comment, "out/embeddings/umap.tsv.gz")
write_tsv(umap_coords, "out/embeddings/umap.tsv.gz", append=TRUE, col_names=TRUE)

          