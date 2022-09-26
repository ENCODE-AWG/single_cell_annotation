Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(remotes)

remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(DoubletFinder)

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

proj <- readRDS(file = input_paths[["project_in"]])

print(proj) ####

# pK identification (no ground-truth)
sweep.list <- paramSweep_v3(proj, PCs = 1:30, sct = TRUE, num.cores = snakemake@threads)
sweep.stats <- summarizeSweep(sweep.list)
bcmvn <- find.pK(sweep.stats)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic doublet proportion estimate
annotations <- proj@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp.poi <- round(params[["doublet_rate"]] * nrow(proj@meta.data)) 
nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

# run DoubletFinder
proj <- doubletFinder_v3(
    seu = proj, 
    PCs = 1:30, 
    pK = optimal.pk,
    nExp = nExp.poi.adj,
    sct = TRUE
)
print(head(proj@meta.data)) ####
doublets_col <- grep("DF.classifications", names(proj@meta.data), value = TRUE)
proj$doublet_classification <- proj[[doublets_col]]
proj_filtered <- subset(proj, subset = doublet_classification == "Singlet")

print(proj)
print(proj_filtered)
saveRDS(proj, file = output_paths[["project_out_all"]])
saveRDS(proj_filtered, file = output_paths[["project_out_filtered"]])

sink(type = "message")
sink()