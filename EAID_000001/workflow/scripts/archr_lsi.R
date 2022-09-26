Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(parallel)

library(ArchR)

# Bug-patched function
.getRowVars <- function(
  ArrowFiles = NULL,
  seqnames = NULL,
  useMatrix = NULL,
  useLog2 = FALSE,
  threads = 1
  ){
  
  .combineVariances <- function(dfMeans = NULL, dfVars = NULL, ns = NULL){
    #https://rdrr.io/cran/fishmethods/src/R/combinevar.R
    if(ncol(dfMeans) != ncol(dfVars) | ncol(dfMeans) != length(ns)){
      stop("Means Variances and Ns lengths not identical")
    }
    #Check if samples have NAs due to N = 1 sample or some other weird thing.
    #Set it to min non NA variance
    dfVars <- lapply(seq_len(nrow(dfVars)), function(x){
      vx <- dfVars[x, , drop = FALSE] # https://github.com/GreenleafLab/ArchR/commit/d28fa1d61fd3a3444a01721785b6f32bce8cbf48
      if(any(is.na(vx))){
        vx[is.na(vx)] <- min(vx[!is.na(vx)])
      }
      vx
    }) %>% Reduce("rbind", .)
    combinedMeans <- rowSums(t(t(dfMeans) * ns)) / sum(ns)
    summedVars <- rowSums(t(t(dfVars) * (ns - 1)) + t(t(dfMeans^2) * ns))
    combinedVars <- (summedVars - sum(ns)*combinedMeans^2)/(sum(ns)-1)
    data.frame(combinedVars = combinedVars, combinedMeans = combinedMeans)
  }
  featureDF <- ArchR:::.getFeatureDF(ArrowFiles, useMatrix)
  if(!is.null(seqnames)){
    featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnames),]
  }
  rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))
  fnames <- rownames(featureDF)
  featureDF <- S4Vectors::split(featureDF, as.character(featureDF$seqnames))
  ns <- lapply(seq_along(ArrowFiles), function(y){
    length(ArchR:::.availableCells(ArrowFiles[y], useMatrix))
  }) %>% unlist
  #Compute RowVars
  summaryDF <- ArchR:::.safelapply(seq_along(featureDF), function(x){
    
    o <- h5closeAll()
    seqx <- names(featureDF)[x]
    meanx <- matrix(NA, ncol = length(ArrowFiles), nrow = nrow(featureDF[[x]]))
    varx <- matrix(NA, ncol = length(ArrowFiles), nrow = nrow(featureDF[[x]]))
    for(y in seq_along(ArrowFiles)){
      if(useLog2){
        meanx[, y] <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqx, "/rowMeansLog2"))
        varx[, y] <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqx, "/rowVarsLog2")) 
      }else{
        meanx[, y] <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqx, "/rowMeans"))
        varx[, y] <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqx, "/rowVars"))
      }
    }
    cbind(featureDF[[x]], DataFrame(.combineVariances(meanx, varx, ns)))
  }, threads = threads) %>% Reduce("rbind", .)
  summaryDF <- summaryDF[fnames, , drop = FALSE]
  
  return(summaryDF)
}
assignInNamespace(".getRowVars", .getRowVars, ns="ArchR")

library(chromVARmotifs)

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

##########

# Conduct ATAC LSI dimensionality reduction
proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix", 
    name = "LSI_ATAC", 
    # iterations = 4, 
    # clusterParams = list( 
    #     resolution = c(0.2,0.2,0.6,0.8), 
    #     sampleCells = 10000, 
    #     n.start = 10
    # ), 
    # varFeatures = 25000,
    # dimsToUse = 1:30,
    logFile = log_paths[["lsi_atac"]]
)

##########

# Calculate UMAP coordinates from ATAC LSI values
proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "LSI_ATAC", 
    name = "UMAP_ATAC", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
)

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