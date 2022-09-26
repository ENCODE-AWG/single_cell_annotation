console_log <- file(snakemake@log[[1]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

in_dir <- snakemake@input[["qc_dir"]]
out_dir <- snakemake@output[["out_dir"]]

dir.create(out_dir)

for (sample in list.files(in_dir)) {
    dir.create(file.path(out_dir, sample))

    in_path_meta <- file.path(in_dir, sample, paste0(sample, "-Pre-Filter-Metadata.rds"))
    print(in_path_meta) ####
    out_path_meta <- file.path(out_dir, sample, "metadata.tsv")

    r <- readRDS(in_path_meta)
    print(r) ####
    write.table(r, out_path_meta, sep = '\t', row.names = FALSE, quote = FALSE)

    in_path_doublets <- file.path(in_dir, sample, paste0(sample, "-Doublet-Summary.rds"))
    out_path_doublets <- file.path(out_dir, sample, "doublets.tsv")

    r <- readRDS(in_path_doublets)
    print(r) ####
    res <- r$doubletResults
    d <- data.frame(
        barcode = row.names(r$originalDataUMAP), 
        doubletScoreUMAP = r$originalDataUMAP$score,
        doubletEnrichUMAP = r$originalDataUMAP$enrichment,
        doubletEnrichLSI = res$doubletEnrichLSI, 
        doubletScoreLSI = res$doubletScoreLSI
    )
    write.table(d, out_path_doublets, sep = '\t', row.names = FALSE, quote = FALSE)
}

sink(type = "message")
sink()