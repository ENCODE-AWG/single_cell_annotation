import gzip

HEADER = """
Cell type labels for each cell
# cell_id: Cell ID used in integrated analysis
# cell_type_id: ENCODE ID of the corresponding cell type object
# cell_type_name: Common name of the cell type
# membership_score: A numeric score for the labeling (Placeholder)
"""

COLUMNS = "cell_id\tcell_type_id\tcell_type_name\tmembership_score\n"

def main(cell_data_path, out_path):
    with open(cell_data_path) as f, gzip.open(out_path, "wt") as fo:
        fo.write(HEADER)
        fo.write(COLUMNS)

        h = f.readline().rstrip('\n').split('\t')
        barcode_ind = 0
        cell_type_ind = h.index("cell_labels")

        for line in f:
            entries = line.rstrip('\n').split('\t')

            barcode = entries[barcode_ind]
            cell_type = entries[cell_type_ind]
            cell_type_id = cell_type.replace(" ", "_")
            
            line = f"{barcode}\t{cell_type_id}\t{cell_type}\tNA\n"
            fo.write(line)

cell_data_path, = snakemake.input

out_path, = snakemake.output

main(cell_data_path, out_path)