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
# head(read.table(file = input_paths[["label_data"]], sep = '\t', header = FALSE)) ####
clustdata <- read.table(file = input_paths[["label_data"]], sep = '\t', header = FALSE, row.names = 1, comment.char = "")
head(clustdata) ####

cellnames <- row.names(clustdata)
proj_index <- match(getCellNames(ArchRProj = proj), cellnames)
proj_index <- proj_index[!is.na(proj_index)]
head(proj_index) ####
length(proj_index) ####
head(clustdata[[1]][proj_index]) ####
proj <- addCellColData(
    ArchRProj = proj,
    data = clustdata[[1]][proj_index],
    name = "seurat_label",
    cells = cellnames[proj_index],
    force = FALSE
)

# Plot ATAC clusters
p1 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "cellColData", 
    name = "seurat_label", 
    embedding = "UMAP_Harmony",
    logFile = log_paths[["umap_plot"]]
)

plotPDF(p1, name = "umap_seurat_label.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

cM <- confusionMatrix(proj$Clusters_ATAC, proj$seurat_label)
cM <- as.matrix(cM)
print(cM) ####
print(colnames(cM)) ####
cM <- cM[, !is.na(colnames(cM))]
cM <- cbind(cM, Unknown = 0.01)
print(cM) ####
labelOld <- rownames(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
proj$cell_labels <- mapLabels(proj$Clusters_ATAC, newLabels = labelNew, oldLabels = labelOld)
unknowns <- proj$cell_labels == "Unknown"
proj$cell_labels[unknowns] <- proj$cell_labels_link[unknowns]

label_data <- getCellColData(ArchRProj = proj)
print(label_data) ####
write.table(label_data, output_paths[["labels"]], quote = FALSE, sep = '\t', col.names = NA)


p <- pheatmap::pheatmap(
    # mat = cM[rowSums(cM)>0,], 
    mat = cM / rowSums(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)
plotPDF(p, name = "seurat_label_cm.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# Plot ATAC clusters
p1 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "cellColData", 
    name = "cell_labels", 
    embedding = "UMAP_Harmony",
    logFile = log_paths[["umap_plot"]]
)

plotPDF(p1, name = "umap_full_label.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

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