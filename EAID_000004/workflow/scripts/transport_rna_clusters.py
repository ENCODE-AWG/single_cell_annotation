
REV_COMP = str.maketrans("ATGC", "TACG")
def reverse_complement(seq):
    return str.translate(seq, REV_COMP)[::-1]

def load_whitelists(wl_atac_path, wl_rna_path):
    bc_map = {}
    with open(wl_atac_path) as a, open(wl_rna_path) as r:
        for line_a, line_r in zip(a, r):
            bc_a = line_a.rstrip("\n")
            bc_r = line_r.rstrip("\n")

            bc_map[bc_r] = bc_a

    return bc_map

def main(in_path, out_path, wl_atac_path, wl_rna_path):
    bc_map = load_whitelists(wl_atac_path, wl_rna_path)

    bcs = []
    with open(in_path) as fi, open(out_path, "w") as fo:
        h = fi.readline().rstrip('\n').split('\t')
        label_ind = h.index("labels_named")

        for line in fi:
            entries = line.rstrip('\n').split('\t')
            bc_in = entries[0]
            label = entries[label_ind]
            dataset, bc_rna = bc_in.rsplit('_', 1)

            if bc_rna in bc_map:
                bc_atac = bc_map[bc_rna]
            else:
                continue

            fo.write(f"{dataset}#{bc_atac}\t{label}\n")
            fo.write(f"{dataset}#{reverse_complement(bc_atac)}\t{label}\n")


in_path = snakemake.input["seurat_data"]
wl_atac_path = snakemake.input["wl_atac"]
wl_rna_path = snakemake.input["wl_rna"]

out_path, = snakemake.output

main(in_path, out_path, wl_atac_path, wl_rna_path)