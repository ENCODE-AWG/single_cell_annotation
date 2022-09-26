import os
import gzip

HEADER = """
# cell_id: Cell ID used in integrated analysis
# rna_dataset: ENCODE snRNA-Seq dataset ID
# rna_barcode: snRNA-Seq barcode
# atac_dataset: ENCODE snATAC-Seq dataset ID
# atac_barcode: snATAC-Seq barcode
# rna_umi_count: snRNA UMI's per cell
# atac_fragment_count: snATAC-Seq fragments per cell
# rna_frac_mito: Fraction of mitochondrial RNA reads
# rna_frac_ribo: Fraction of ribosomal RNA reads
# atac_tss_enrichment: snATAC-Seq transcription start site enrichment
# passed_filtering: A binary (0/1) value indicating whether the cell passed manual filtering
"""

COLUMNS = (
    "cell_id\trna_dataset\trna_barcode\tatac_dataset\tatac_barcode\t"
    "rna_umi_count\tatac_fragment_count\trna_frac_mito\trna_frac_ribo\tatac_tss_enrichment\tpassed_filtering\n"
)

def load_records(metadata_path):
    records = {}
    with open(metadata_path) as f:
        h = f.readline().rstrip('\n').split('\t')
        cell_id_ind = h.index("cellNames")
        frag_ind = h.index("nFrags")
        tss_ind = h.index("TSSEnrichment")

        for line in f:
            entries = line.rstrip('\n').split('\t')

            cell_id = entries[cell_id_ind]
            dataset, barcode = cell_id.split("#")
            dataset_parsed = dataset.split("_")[-1].split("-")[0]
            
            frag_count = int(entries[frag_ind])
            tss_enr = float(entries[tss_ind])

            record = {
                "dataset": dataset_parsed,
                "barcode": barcode,
                "frag_count": frag_count,
                "tss_enr": tss_enr,
                "pass_filter": False
            }
            records[cell_id] = record

    return records

def load_final_data(final_data_path):
    ids = []
    with open(final_data_path) as f:
        next(f)
        for line in f:
            cell_id = line.rstrip('\n').split('\t')[0]
            ids.append(cell_id)

    return ids

def main(metadata_dir, final_data_path, out_path):
    metadata_paths = [os.path.join(metadata_dir, i, "metadata.tsv") for i in os.listdir(metadata_dir) if not i.startswith(".")]
    records = {}
    for i in metadata_paths:
        records |= load_records(i)

    final_ids = load_final_data(final_data_path)
    for i in final_ids:
        records[i]["pass_filter"] = True

    with gzip.open(out_path, "wt") as f:
        f.write(HEADER)
        f.write(COLUMNS)

        for cell_id, r in records.items():
            line = f"{cell_id}\tNA\tNA\t{r['dataset']}\t{r['barcode']}\tNA\t{r['frag_count']}\tNA\tNA\t{r['tss_enr']}\t{int(r['pass_filter'])}\n"
            f.write(line)

metadata_dir = snakemake.input["metadata"]
final_data_path = snakemake.input["final_data"]

out_path, = snakemake.output

main(metadata_dir, final_data_path, out_path)