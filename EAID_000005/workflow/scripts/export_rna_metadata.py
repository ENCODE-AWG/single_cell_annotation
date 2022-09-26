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
        barcode_ind = 0
        dataset_ind = h.index("orig.ident")
        count_ind = h.index("nCount_RNA")
        mito_ind = h.index("percent.mt")
        ribo_ind = h.index("percent.ribo")

        for line in f:
            entries = line.rstrip('\n').split('\t')

            dataset = entries[dataset_ind]
            barcode = entries[barcode_ind]
            cell_id = f"{dataset}_{barcode}"
            dataset_parsed = dataset.split("_")[0].split("-")[0]
            
            umi_count = int(entries[count_ind])
            frac_mito = float(entries[mito_ind]) / 100
            frac_ribo = float(entries[ribo_ind]) / 100

            record = {
                "dataset": dataset_parsed,
                "barcode": barcode,
                "umi_count": umi_count,
                "frac_mito": frac_mito,
                "frac_ribo": frac_ribo,
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

def main(metadata_paths, final_data_path, out_path):
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
            line = f"{cell_id}\t{r['dataset']}\t{r['barcode']}\tNA\tNA\t{r['umi_count']}\tNA\t{r['frac_mito']}\t{r['frac_ribo']}\tNA\t{int(r['pass_filter'])}\n"
            f.write(line)

metadata_paths = snakemake.input["metadata"]
final_data_path = snakemake.input["final_data"]

out_path, = snakemake.output

main(metadata_paths, final_data_path, out_path)