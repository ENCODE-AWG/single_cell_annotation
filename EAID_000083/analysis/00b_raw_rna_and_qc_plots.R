suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(SeuratDisk)
  library(patchwork)
  library(data.table)
  library(GenomicRanges)
})
# This is adapted from 05_analysis/meetings/Mar30/PANC_analysis.R
# The aim is just to produce a list of Seurat objects with the filtered cells + 
# raw counts matrices (subsetting to protein_coding, lincrna, and TR/IG genes)
source("analysis/utils.R")

meta <- read_tsv("config/panc_samples.tsv")
cutoffs <-read_tsv("config/sample_cutoffs.tsv")

panc_all <- meta %>% pull(id)
panc_multiome <- meta %>% filter(multiome_series != "n/a") %>% pull(id)


qc <- fread("out/metadata.tsv.gz", skip="cell_id\tsample_name") 


tss_plots <- map(panc_all, function(s) {
  c <- filter(cutoffs, sample_id == s)
  cutoff_plot_hexbin(filter(qc, sample_name == s), log10(atac_fragment_count), atac_tss_enrichment, log10(c$min_nFrags), 5) +
    theme_bw() + ggtitle(s) + ylim(0, 30)
})
tss_plot <- patchwork::wrap_plots(tss_plots, ncol=3)

rna_count_plots <- map(panc_all, function(s) {
  c <- filter(cutoffs, sample_id == s)
  read_count_cutoff_plot(filter(qc, sample_name == s) %>% pull(rna_umi_count), c$min_umis) +
    theme_bw() + ggtitle(s)
})
rna_count_plot <- patchwork::wrap_plots(rna_count_plots, ncol=3)

gene_metadata <- read_tsv("data/ENCODE/samples/W61_PANC/rna/features.tsv.gz", col_names = c("id", "symbol", "type"))

gencode <- read_gtf("data/gencode.v29.annotation.gtf.gz", attributes=c("gene_type", "gene_id", "gene_name"), skip="chr1") %>%
  as_tibble() %>%
  filter(feature == "gene")

keeper_genes <- gencode %>%
  filter(gene_type %in% c("protein_coding", "lincRNA") | str_detect(gene_type, "^IG_|TR_"))


rna_mats <- map(panc_all, function(s) {
  c <- filter(cutoffs, sample_id == s)
  passing_cells <- rna_qc %>%
    filter(umis >= c$min_umis, sample_id == s) %>%
    pull(cell_barcode)
  passing_cells_atac <- atac_qc %>% 
    filter(TSSEnrichment >= 5, sample_id == s, nFrags >= c$min_nFrags) %>%
    pull(cell_barcode)
  if (filter(meta, id == s)[["multiome_series"]] != "n/a")
    passing_cells <- intersect(passing_cells, passing_cells_atac)
  mat <- readRDS(sprintf("data/seurat/%s/rna_raw.rds", s))[keeper_genes$gene_id,passing_cells]
  rownames(mat) <- keeper_genes$gene_name
  mat
})

seurat_raw <- mclapply(seq_along(rna_mats), mc.cores=4, function(i) {
  mat <- rna_mats[[i]]
  sample <- panc_all[i]
  meta.data <- data.frame(sample=rep_len(sample, ncol(mat)))
  rownames(meta.data) <- colnames(mat)
  CreateSeuratObject(mat, project=sample, meta.data=meta.data)
})

saveRDS(seurat_raw, "data/seurat_raw.rds")

