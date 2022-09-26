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

proj <- readRDS(file = input_paths[["project_in"]])
print(proj) ####

umap_coords <- Embeddings(proj, reduction = "umap")
head(umap_coords) ####
write.table(umap_coords, file=output_paths[["umap_coords"]], quote=FALSE, sep='\t', col.names = NA)

emb_coords <- Embeddings(proj, reduction = "harmony")
head(emb_coords) ####
write.table(emb_coords, file=output_paths[["emb_coords"]], quote=FALSE, sep='\t', col.names = NA)

sink(type = "message")
sink()