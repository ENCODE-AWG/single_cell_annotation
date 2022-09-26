Sys.setenv(CONDA_BUILD_SYSROOT="/")

console_log <- file(snakemake@log[["console"]], open = "wt")
sink(console_log)
sink(console_log, type = "message")

library(Seurat)

params = snakemake@params 
input_paths = snakemake@input
output_paths = snakemake@output
log_paths = snakemake@log

set.seed(params[["seed"]])

dir.create(output_paths[["markers"]])

proj <- readRDS(file = input_paths[["project_in"]])

label_names <- unique(proj$labels_named)

for (i in label_names){
    markers <- FindMarkers(proj, ident.1 = i, group.by = 'labels_named')
    out_path <- file.path(output_paths[["markers"]], paste0(gsub(" ", "_", i), ".tsv"))
    write.table(markers, file=out_path, quote=FALSE, sep='\t', col.names = NA)
}

sink(type = "message")
sink()