import os
import gzip

HEADER = """
# Marker genes for {label_name}
gene_id: ENSEMBL Gene ID
gene_name: Common name of the gene
avg_log2FC: Log2 fold-chage of the average expression between the given cell type and the remaining cells (min 0.25)
p_val: Raw p-value using Wilcoxon Rank Sum test 
p_val_adj: Bonferroni-adjusted p-value across genes
is_enriched: A binary (0/1) value indicating whether the given gene is enriched for the cell type (adjusted p < 0.05)
"""

COLUMNS = "gene_id\tgene_name\tavg_log2FC\tavg_log2FC\tp_val\tp_val_adj\tis_enriched\n"

def build_gene_map(genes_paths):
    gene_map = {}
    for p in genes_paths:
        with open(p) as f:
            for line in f:
                gene_id, gene_name = line.rstrip('\n').split('\t')[:2]
                gene_map[gene_name] = gene_id

    return gene_map

def export_label(in_path, out_path, gene_map, label_name):
    with open(in_path) as f, gzip.open(out_path, "wt") as fo:
        fo.write(HEADER.format(label_name=label_name))
        fo.write(COLUMNS)

        h = f.readline().rstrip('\n').split('\t')
        gene_name_ind = 0
        lfc_ind = h.index("avg_log2FC")
        pval_ind = h.index("p_val")
        pval_adj_ind = h.index("p_val_adj")

        for line in f:
            entries = line.rstrip('\n').split('\t')

            gene_name = entries[gene_name_ind]
            avg_logFC = float(entries[lfc_ind])
            p_val = float(entries[pval_ind])
            p_val_adj = float(entries[pval_adj_ind])

            gene_id = gene_map[gene_name]

            is_enriched = int((p_val_adj < 0.05) and (avg_logFC > 0))
            
            line = f"{gene_id}\t{gene_name}\t{avg_logFC}\t{p_val}\t{p_val_adj}\t{is_enriched}\n"
            fo.write(line)

def main(markers_dir, genes_paths, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    gene_map = build_gene_map(genes_paths)
    labels = [i for i in os.listdir(markers_dir) if not i.startswith(".")]
    for l in labels:
        markers_path = os.path.join(markers_dir, l)
        label_id = os.path.splitext(l)[0]
        label_name = label_id.replace("_", " ")
        out_path = os.path.join(out_dir, f"{label_id}.tsv.gz")
        export_label(markers_path, out_path, gene_map, label_name)

markers_dir = snakemake.input["markers"]
genes_paths = snakemake.input["genes"]

out_dir, = snakemake.output

main(markers_dir, genes_paths, out_dir)