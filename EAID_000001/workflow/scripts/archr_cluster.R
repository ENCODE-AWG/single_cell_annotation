Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(parallel)

library(ArchR)
library(chromVARmotifs)
library(Seurat)

addArchRVerbose(verbose = FALSE)
addArchRChrPrefix(chrPrefix = FALSE)

# Disable HDF5 file locking
# Workaround for HDF5 I/O issues on NFS
# https://github.com/GreenleafLab/ArchR/issues/248#issuecomment-789453997
Sys.setenv("HDF5_USE_FILE_LOCKING" = "FALSE")
Sys.setenv("RHDF5_USE_FILE_LOCKING" = "FALSE")

##########

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
threads = snakemake@threads
log_paths = snakemake@log

seed <- params[["seed"]]
set.seed(seed)

addArchRThreads(threads = threads)

##########

# Load and move project
proj <- loadArchRProject(path = input_paths[["project_in"]], force = FALSE, showLogo = FALSE)
proj <- saveArchRProject(
    ArchRProj = proj,
    outputDirectory = output_paths[["project_out"]],
    overwrite = TRUE,
    load = TRUE,
    logFile = log_paths[["move"]],
)

# Cluster cells by ATAC LSI values
proj <- addClusters(
    input = proj,
    reducedDims = "Harmony_ATAC",
    method = "Seurat",
    name = "Clusters_ATAC",
    resolution = 0.8,
    logFile = log_paths[["cluster_atac"]]
)

##########

# Plot ATAC clusters
p1 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "cellColData", 
    name = "Clusters_ATAC", 
    embedding = "UMAP_Harmony",
    logFile = log_paths[["umap_plot"]]
)

plotPDF(p1, name = "umap_atac_clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

##########

saveArchRProject(
    ArchRProj = proj,
    outputDirectory = output_paths[["project_out"]],
    overwrite = TRUE,
    load = FALSE,
    logFile = log_paths[["save"]],
)


sink(type = "message")
sink()