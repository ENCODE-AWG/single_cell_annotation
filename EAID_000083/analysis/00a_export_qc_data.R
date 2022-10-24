suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(dtplyr)
  library(parallel)
  library(Matrix)
})
source("analysis/utils.R")

meta <- read_tsv("config/panc_samples.tsv", col_types="cccccc")
cutoffs <-read_tsv("config/sample_cutoffs.tsv")

panc_all <- meta$id
panc_multiome <- meta %>% filter(multiome_series != "n/a") %>% pull(id)

gencode <- read_gtf("data/gencode.v29.annotation.gtf.gz", attributes=c("gene_type", "gene_id", "gene_name"), skip="chr1") %>%
  as_tibble() %>%
  filter(feature == "gene")
mito_genes <- gencode %>% filter(seqnames == "chrM") %>% pull(gene_id)
ribo_genes <- gencode %>% filter(gene_type == "rRNA") %>% pull(gene_id)

# rna_dataset: ENCODE scRNA-Seq dataset ID
# rna_barcode: scRNA-Seq barcode
# rna_frac_mito: Fraction of mitochondrial RNA reads
# rna_frac_ribo: Fraction of ribosomal RNA reads
# rna_umi_count: scRNA UMI’s per cell
rna_qc <- mclapply(panc_all, mc.cores=6, function(s) {
  mat <- readRDS(sprintf("../../04_data/seurat/%s/rna_raw.rds", s))
  tibble(
    sample_name = s,
    rna_barcode = colnames(mat),
    cell_id = sprintf("%s#%s", sample_name, rna_barcode),
    rna_umi_count = colSums(mat),
    rna_frac_mito = colSums(mat[mito_genes,]) / rna_umi_count,
    rna_frac_ribo = colSums(mat[ribo_genes,]) / rna_umi_count
  )
}) %>% bind_rows()

# rna_dataset: ENCODE scRNA-Seq dataset ID
rna_qc <- rna_qc %>%
  inner_join(
    select(meta, rna_dataset=rna_id, sample_name=id, multiome_series),
    by = "sample_name"
  ) %>%
  mutate(multiome_series=na_if(multiome_series, "n/a"))

# atac_fragment_count: snATAC-Seq fragments per cell
# atac_tss_enrichment: snATAC-Seq transcription start site enrichment
atac_qc <- mclapply(panc_all, mc.cores=6, function(s) {
  atac_qc <- readRDS(sprintf("data/archr/%1$s/ArchR_Project_unfiltered/QualityControl/%1$s/%1$s-Pre-Filter-Metadata.rds", s))
  tibble(
    sample_name = s,
    cell_id = atac_qc$cellNames,
    atac_fragment_count = atac_qc$nFrags,
    atac_tss_enrichment = atac_qc$TSSEnrichment
  )
}) %>% bind_rows()

# atac_barcode: snATAC-Seq barcode
atac_qc <- atac_qc %>%
  mutate(rna = str_remove(cell_id, "^[^#]+#")) %>%
  left_join(read_tsv("config/10x_barcode_translation.txt.gz", col_types="cc"), by="rna") %>%
  mutate(atac_barcode = if_else(is.na(atac), rna, atac)) %>%
  select(!c(rna, atac))
# atac_dataset: ENCODE snATAC-Seq dataset ID
atac_qc <- atac_qc %>%
  inner_join(
    select(meta, atac_dataset=atac_id, sample_name=id, multiome_series),
    by = "sample_name"
  ) %>%
  mutate(multiome_series=na_if(multiome_series, "n/a"))

joint_qc <- mclapply(panc_all, mc.cores=6, function(s) {
  m <- meta %>% filter(id == s)
  rna_qc <- filter(rna_qc, rna_dataset == m$rna_id)
  atac_qc <- filter(atac_qc, atac_dataset == m$atac_id)
  if (m$multiome_series == "n/a") {
    bind_rows(rna_qc, atac_qc)
  } else {
    rna_qc %>%
      lazy_dt() %>%
      select(!multiome_series) %>%
      full_join(atac_qc, by=c("sample_name", "cell_id")) %>%
      select(!multiome_series) %>%
      mutate(
        across(c(rna_umi_count, rna_frac_mito, rna_frac_ribo, atac_fragment_count, atac_tss_enrichment),
               ~ replace_na(.x, 0))
      ) %>%
      mutate(
        rna_dataset = m$rna_id,
        atac_dataset = m$atac_id
      ) %>%
      collect()
  }
}) %>% bind_rows()


final_qc <- joint_qc %>%
  inner_join(cutoffs, by=c("sample_name"="sample_id")) %>%
  mutate(
    passed_atac = atac_fragment_count >= min_nFrags & atac_tss_enrichment >= 5,
    passed_rna = rna_umi_count >= min_umis,
    passed_filtering = (passed_atac | is.na(atac_dataset)) & (passed_rna | is.na(rna_dataset))
  ) %>% 
  select(cell_id, sample_name, rna_dataset, rna_barcode, atac_dataset, atac_barcode, 
         rna_umi_count, atac_fragment_count, 
         rna_frac_mito, rna_frac_ribo, atac_tss_enrichment, 
         passed_filtering)

header_text <- c(
  "# cell_id: Cell ID used in integrated analysis",
  "# rna_dataset: ENCODE snRNA-Seq dataset ID",
  "# rna_barcode: snRNA-Seq barcode",
  "# atac_dataset: ENCODE snATAC-Seq dataset ID",
  "# atac_barcode: snATAC-Seq barcode",
  "# rna_umi_count: snRNA UMI's per cell",
  "# atac_fragment_count: snATAC-Seq fragments per cell",
  "# rna_frac_mito: Fraction of mitochondrial RNA reads",
  "# rna_frac_ribo: Fraction of ribosomal RNA reads",
  "# atac_tss_enrichment: snATAC-Seq transcription start site enrichment",
  "# passed_filtering: A binary (0/1) value indicating whether the cell passed manual filtering"
)


write_lines(header_text, "out/metadata.tsv.gz")
write_tsv(final_qc, "out/metadata.tsv.gz", append=TRUE, col_names=TRUE)
