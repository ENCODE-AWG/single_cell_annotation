import gzip

HEADER_EMB = """
# Unified cell embeddings integrated using Harmony
# cell_id: Cell ID used in integrated analysis
# harmony_1, harmony_2, … harmony_50: Columns of a 50-dimensional Harmony embedding vector
"""

HEADER_UMAP = """
# UMAP coordinates of each cell
# cell_id: Cell ID used in integrated analysis
# UMAP_1, UMAP_2: UMAP x and y coordinates, respectively
"""

def main(emb_in_path, umap_in_path, emb_out_path, umap_out_path):
    with open(emb_in_path) as f, gzip.open(emb_out_path, "wt") as fo:
        fo.write(HEADER_EMB)

        h = f.readline()
        fo.write("cell_id" + h)

        for line in f:
            fo.write(line)

    with open(umap_in_path) as f, gzip.open(umap_out_path, "wt") as fo:
        fo.write(HEADER_UMAP)

        h = f.readline()
        fo.write("cell_id" + h)

        for line in f:
            fo.write(line)

emb_in_path = snakemake.input["emb"]
umap_in_path = snakemake.input["umap"]

emb_out_path = snakemake.output["emb"]
umap_out_path = snakemake.output["umap"]

main(emb_in_path, umap_in_path, emb_out_path, umap_out_path)