Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(parallel)
library(BiocGenerics)

library(ArchR)
library(chromVARmotifs)

addArchRVerbose(verbose = FALSE)
addArchRChrPrefix(chrPrefix = FALSE)

# Disable HDF5 file locking
# Workaround for HDF5 I/O issues on NFS
# https://github.com/GreenleafLab/ArchR/issues/248#issuecomment-789453997
Sys.setenv("HDF5_USE_FILE_LOCKING" = "FALSE")
Sys.setenv("RHDF5_USE_FILE_LOCKING" = "FALSE")

partial <- function(f, ...) {
    l <- list(...)
    function(...) {
        do.call(f, c(l, list(...)))
    }
}

plot_sample <- function(proj, log_path, sample_name) {
    idxSample <- BiocGenerics::which(proj$Sample %in% sample_name)
    cellsSample <- proj$cellNames[idxSample]
    proj <- proj[cellsSample, ]
    p <- plotEmbedding(
        ArchRProj = proj, 
        colorBy = "cellColData", 
        name = "Clusters_WNN", 
        embedding = "UMAP_Harmony",
        logFile = log_path
    )
    plt_name <- sprintf("umap_clusters_harmony_%s.pdf", sample_name)
    plotPDF(p, name = plt_name, ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
    p
}

main <- function(params, input_paths, output_paths, threads, log_paths) {
    seed <- params[["seed"]]
    set.seed(seed)

    addArchRThreads(threads = threads)

    # Load and move project
    proj <- loadArchRProject(path = input_paths[["project_in"]], force = FALSE, showLogo = FALSE)
    proj <- saveArchRProject(
        ArchRProj = proj,
        outputDirectory = output_paths[["project_out"]],
        overwrite = TRUE,
        load = TRUE,
        logFile = log_paths[["move"]],
    )

    p <- plotEmbedding(
        ArchRProj = proj, 
        colorBy = "cellColData", 
        name = "Sample", 
        embedding = "UMAP_ATAC",
        logFile = log_paths[["umap_plot"]]
    )
    plotPDF(p, name = "umap_pre_harmony_datasets.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

    # Conduct Harmony batch correction
    proj <- addHarmony(
        ArchRProj = proj,
        reducedDims = "LSI_ATAC",
        name = "Harmony_ATAC",
        groupBy = "Sample",
        verbose = FALSE
    )

    # Calculate UMAP coordinates from Harmony-adjusted values
    proj <- addUMAP(
        ArchRProj = proj, 
        reducedDims = "Harmony_ATAC", 
        name = "UMAP_Harmony", 
        nNeighbors = 30, 
        minDist = 0.5, 
        metric = "cosine",
    )

    # plot_fn <- partial(plot_sample, proj = proj, log_path = log_paths[["umap_plot"]])
    # plts <- lapply(params[["sample_names"]], plot_fn)

    p <- plotEmbedding(
        ArchRProj = proj, 
        colorBy = "cellColData", 
        name = "Sample", 
        embedding = "UMAP_Harmony",
        logFile = log_paths[["umap_plot"]]
    )
    plotPDF(p, name = "umap_harmony_datasets.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
    


    # p2 <- plotEmbedding(
    #     ArchRProj = proj, 
    #     colorBy = "cellColData", 
    #     name = "ClustersHarmony", 
    #     embedding = "UMAPHarmony",
    #     logFile = log_paths[["umap_plot"]]
    # )
    # ggAlignPlots(p1, p2, type = "h")
    # plotPDF(p1, p2, name = "umap_clusters_harmony.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

    saveArchRProject(
        ArchRProj = proj,
        outputDirectory = output_paths[["project_out"]],
        overwrite = TRUE,
        load = FALSE,
        logFile = log_paths[["save"]],
    )

}

main(snakemake@params, snakemake@input, snakemake@output, snakemake@threads, snakemake@log)

sink(type = "message")
sink()