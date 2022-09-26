Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(parallel)

library(ArchR)
library(pheatmap)

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

# Load project
proj <- loadArchRProject(path = input_paths[["project_in"]], force = FALSE, showLogo = FALSE)

umap_coords <- getEmbedding(ArchRProj = proj, embedding = "UMAP_Harmony", returnDF = TRUE)
write.table(umap_coords, output_paths[["umap_coords"]], quote = FALSE, sep = '\t', col.names = NA)

emb_coords <- getReducedDims(
  ArchRProj = proj,
  reducedDims = "Harmony_ATAC",
  returnMatrix = TRUE,
  dimsToUse = NULL,
  scaleDims = NULL,
  corCutOff = 0.75
)
write.table(emb_coords, output_paths[["emb_coords"]], quote = FALSE, sep = '\t', col.names = NA)

markersGS <- getMarkerFeatures(
  ArchRProj = proj,
  groupBy = "cell_labels",
  useGroups = NULL,
  bgdGroups = NULL,
  useMatrix = "GeneScoreMatrix",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  normBy = NULL,
  testMethod = "wilcoxon",
  maxCells = 500,
  scaleTo = 10^4,
  threads = getArchRThreads(),
  k = 100,
  bufferRatio = 0.8,
  binarize = FALSE,
  useSeqnames = NULL,
  verbose = TRUE,
  logFile = log_paths[["markers"]]
)

markerList <- getMarkers(
  seMarker = markersGS,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  n = NULL,
  returnGR = FALSE
)

dir.create(output_paths[["markers"]])
for (label in names(markerList)) {
    markers <- markerList[[label]]
    out_path <- file.path(output_paths[["markers"]], paste0(gsub(" ", "_", label), ".tsv"))
    write.table(markers, file=out_path, quote=FALSE, sep='\t', col.names = NA)
}

sink(type = "message")
sink()