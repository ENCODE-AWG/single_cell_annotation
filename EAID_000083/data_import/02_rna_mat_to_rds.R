# wd <- getwd()
# setwd("/oak/stanford/groups/wjg/bparks/multiome/")
# source("renv/activate.R")
# setwd(wd)

args <- commandArgs(trailing=TRUE)
stopifnot(length(args) == 2)
input_path <- args[1]
output_path <- args[2]

suppressPackageStartupMessages({
    library(Seurat)
})

x <- Seurat::Read10X(input_path, gene.column=1)
saveRDS(x, output_path)