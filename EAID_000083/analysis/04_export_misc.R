suppressPackageStartupMessages({
  library(tidyverse)
})

meta <- read_tsv("config/panc_samples.tsv", col_types="cccccc")

datasets <- c(meta$rna_id, meta$atac_id[meta$multiome_series != "n/a"])

write_lines(datasets, "out/datasets.txt")

# Output marker genes

azimuth_markers <- read_tsv("config/azimuth_panc_markers.tsv", col_types="")


# My table columns:
# - gene_id: ENSEMBL Gene ID
# - gene_name: Common name of the gene
# - is_enriched: A binary (0/1) value indicating whether the given gene is enriched for the cell type
# - source: (hubmap | azimuth | hubmap+azimuth)

# From asct+b
hubmap_markers <- list(
  "acinar" = c("AMY2A", "PRSS1", "CTRB1", "CELA3A", "PNLIP", "KRT18", "KRT8"), 
  "ductal" = c("CFTR", "KRT19", "KRT7"), 
  "beta" = c("CHGA", "SYP", "UCHL1", "INS", "ENTPD3"), 
  "alpha_gamma_epsilon" = c("CHGA", "SYP", "GCG", "ARX", "SST"), 
  "t" = c("PTPRC", "CD3E", "CD8A", "CD4"),
  "myeloid" = c("PTPRC", "CD163", "MRC1"),
  "immune" = c("PTPRC", "ITGAX"),
  "stellate" = c("ACTA2", "LAMA2", "PDGFRB")
)

azimuth_markers <- read_tsv("config/azimuth_panc_markers.tsv", col_types="")
azimuth_markers_list <- str_split(azimuth_markers$Markers, ", ") %>% 
  set_names(azimuth_markers$Label)

azimuth_markers_clean <- list(
  "acinar" = azimuth_markers_list$acinar,
  "ductal" = azimuth_markers_list$ductal,
  "beta" = azimuth_markers_list$beta,
  "alpha_gamma_epsilon" = c(azimuth_markers_list$alpha, azimuth_markers_list$gamma, azimuth_markers_list$epsilon),
  "gamma" = azimuth_markers_list$gamma,
  "immune" = azimuth_markers_list$immune,
  "stellate" = c(azimuth_markers_list$activated_stellate, azimuth_markers_list$quiescent_stellate),
  "schwann" = azimuth_markers_list$schwann,
  "endothelial" = azimuth_markers_list$endothelial
)

gene_metadata <- read_tsv("data/ENCODE/samples/W61_PANC/rna/features.tsv.gz", col_names = c("gene_id", "gene_name", "type"))


marker_list <- bind_rows(
  enframe(hubmap_markers, name="cell_type", value="gene_name") %>% mutate("source"= "hubmap"),
  enframe(azimuth_markers_clean, name="cell_type", value="gene_name") %>% mutate("source"="azimuth")
) %>%
  unnest_longer(gene_name) %>%
  inner_join(select(gene_metadata, gene_id, gene_name), by="gene_name") %>%
  group_by(cell_type, gene_name, gene_id) %>%
  summarize(source=paste(sort(source), collapse="+"), .groups="drop")

# Actually save the marker gene files
dir.create("out/markers")
comments <- c(
  "# gene_id: ENSEMBL Gene ID",
  "# gene_name: Common name of the gene",
  "# is_enriched: A binary (0/1) value indicating whether the given gene is enriched for the cell type (FDR < 0.1)",
  "# source: hubmap asct+b (https://doi.org/10.48539/HBM557.XHHV.557) or azimuth (https://azimuth.hubmapconsortium.org/references/#Human%20-%20Pancreas)"
)
for (type in unique(marker_list$cell_type)) {
  file_path <- sprintf("out/markers/%s.tsv.gz", type)
  write_lines(comments, file_path)
  marker_list %>%
    filter(cell_type == type) %>%
    mutate(is_enriched = 1) %>%
    select(gene_id, gene_name, is_enriched, source) %>%
    write_tsv(file_path, append=TRUE, col_names=TRUE)
}

  
