args <- commandArgs(trailing=TRUE)
stopifnot(length(args) == 2)
input_h5ad <- args[1]
output_h5seurat <- args[2]

# source("renv/activate.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SeuratDisk)
})

Convert(input_h5ad, output_h5seurat, 
        overwrite = TRUE,
        verbose = FALSE)

proj <- LoadH5Seurat(
  output_h5seurat,
  assays = c("RNA"),
  misc = FALSE,
  tools = FALSE
)
