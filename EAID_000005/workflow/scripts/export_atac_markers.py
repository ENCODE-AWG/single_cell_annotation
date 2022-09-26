import os
import gzip

HEADER = """
# Marker genes for {label_name}
gene_id: ENSEMBL Gene ID
gene_name: Common name of the gene
avg_log2FC: Log2 fold-chage of the average expression between the given cell type and the remaining cells (min 0.5)
FDR: p-value-derived FDR value for gene 
is_enriched: A binary (0/1) value indicating whether the given gene is enriched for the cell type (FDR < 0.1)
"""

COLUMNS = "gene_id\tgene_name\tavg_log2FC\tFDR\tis_enriched\n"

def load_gtf(gtf_path):
    gtf_data = {}
    with gzip.open(gtf_path, "rt") as i:
        for line in i:
            if line.startswith("#"):
                continue

            chrom, source, feature, start, end, score, strand, frame, attributes = line.rstrip('\n').split('\t')

            if feature != "gene":
                continue

            attributes = [i.strip().split(" ") for i in attributes.rstrip(";").split(";")]
            attributes = {k: v.strip("\"") for k, v in attributes}            
            gene_id = attributes["gene_id"]
            gene_name = attributes.get("gene_name", gene_id)

            gtf_data[gene_name] = gene_id

    return gtf_data


def export_label(in_path, out_path, gene_data, label_name):
    with open(in_path) as f, gzip.open(out_path, "wt") as fo:
        fo.write(HEADER.format(label_name=label_name))
        fo.write(COLUMNS)

        h = f.readline().rstrip('\n').split('\t')
        gene_name_ind = h.index("name")
        lfc_ind = h.index("Log2FC")
        fdr_ind = h.index("FDR")

        for line in f:
            entries = line.rstrip('\n').split('\t')

            gene_name = entries[gene_name_ind]
            avg_logFC = float(entries[lfc_ind])
            fdr = float(entries[fdr_ind])

            gene_id = gene_data.get(gene_name, gene_name)

            is_enriched = int((fdr < 0.1) and (avg_logFC > 0))
            
            line = f"{gene_id}\t{gene_name}\t{avg_logFC}\t{fdr}\t{is_enriched}\n"
            fo.write(line)

def main(markers_dir, out_dir, gtf_path):
    os.makedirs(out_dir, exist_ok=True)
    labels = [i for i in os.listdir(markers_dir) if not i.startswith(".")]
    gene_data = load_gtf(gtf_path)
    for l in labels:
        markers_path = os.path.join(markers_dir, l)
        label_id = os.path.splitext(l)[0]
        label_name = label_id.replace("_", " ")
        out_path = os.path.join(out_dir, f"{label_id}.tsv.gz")
        export_label(markers_path, out_path, gene_data, label_name)

markers_dir = snakemake.input["markers"]
gtf_path = snakemake.input["gtf"]

out_dir, = snakemake.output

main(markers_dir, out_dir, gtf_path)