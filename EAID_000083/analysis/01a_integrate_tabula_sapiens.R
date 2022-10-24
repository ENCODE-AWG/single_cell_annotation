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

seurat_vanilla <- readRDS("data/seurat_raw.rds")

ts_panc <-LoadH5Seurat(
  "data/TS_Pancreas.h5seurat",
  assays = c("RNA"),
  misc = FALSE,
  tools = FALSE
)
ts_panc <- ts_panc %>%
  NormalizeData() %>%
  FindVariableFeatures()

# About 80 seconds
system.time(
  anchors <- mclapply(seurat_vanilla, mc.cores = 8, function(p) {
    features <- SelectIntegrationFeatures(list(ts_panc, p))
    anchors <- FindTransferAnchors(ts_panc, p, features = features)
    anchors
  })
)

cell_labels <- map(seq_along(anchors), function(i) tibble(
    cell_label=TransferData(anchors[[i]], ts_panc$cell_ontology_class)$predicted.id,
    sample = seurat_vanilla[[i]]$sample,
    cell_id = colnames(seurat_vanilla[[i]])
  )) %>%
  bind_rows() 

write_tsv(cell_labels, "data/tabula_sapiens_cell_labels.tsv")
