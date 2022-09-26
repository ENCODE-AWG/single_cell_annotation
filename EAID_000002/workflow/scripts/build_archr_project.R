Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(parallel)

library(ArchR)

addArchRVerbose(verbose = FALSE)
addArchRChrPrefix(chrPrefix = FALSE)

# Disable HDF5 file locking
# Workaround for HDF5 I/O issues on NFS
# https://github.com/GreenleafLab/ArchR/issues/248#issuecomment-789453997
Sys.setenv("HDF5_USE_FILE_LOCKING" = "FALSE")
Sys.setenv("RHDF5_USE_FILE_LOCKING" = "FALSE")

build_archr_project <- function(params, input_paths, output_paths, threads, log_paths) {
    arrow_sample_names <- unlist(params[["sample_names"]])
    seed <- params[["seed"]]
    min_frags <- params[["min_frags"]]
    min_tss_enr <- params[["min_tss_enr"]]

    set.seed(seed)
    
    addArchRThreads(threads = threads)

    addArchRGenome("hg38")

    frag_paths <- unlist(input_paths[["frags"]])

    arrow_output_dir <- output_paths[["arrow_dir"]]
    arrow_output_names <- file.path(arrow_output_dir, arrow_sample_names)
    # print(arrow_output_dir) ####
    dir.create(arrow_output_dir, recursive = TRUE)
    arrows <- createArrowFiles(
        inputFiles = frag_paths,
        sampleNames = arrow_sample_names,
        outputNames = arrow_output_names,
        offsetPlus = 0,
        offsetMinus = 0,
        minFrags = min_frags,
        minTSS = min_tss_enr,
        addTileMat = TRUE,
        addGeneScoreMat = TRUE,
        force = TRUE,
        subThreading = FALSE, # required for no file locking
        logFile = log_paths[["arrow_create"]],
        QCDir = output_paths[["qc_dir"]]
    )

    # Calculate doublet scores
    doub_scores <- addDoubletScores(
        input = arrows,
        k = 10, 
        knnMethod = "UMAP", 
        LSIMethod = 1,
        outDir = output_paths[["qc_dir"]],
        logFile = log_paths[["doublets"]],
    )

    # Create project
    dir.create(output_paths[["project_dir"]])
    proj <- ArchRProject(
        ArrowFiles = arrows, 
        outputDirectory = output_paths[["project_dir"]],
        copyArrows = FALSE,
    )
    
    # Filter doublets
    proj <- filterDoublets(proj)

    saveArchRProject(
        ArchRProj = proj,
        overwrite = TRUE,
        load = FALSE,
        logFile = log_paths[["save"]],
    )

}

build_archr_project(snakemake@params, snakemake@input, snakemake@output, snakemake@threads, snakemake@log)

sink(type = "message")
sink()