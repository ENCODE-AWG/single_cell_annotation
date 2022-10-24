import sys
import gzip

# Usage: python translate_barcodes.py < fragments.tsv.gz | bgzip -c > renamed.tsv.gz
# then run tabix on the bgzipped output

atac_to_rna = {}
for l in gzip.open("config/10x_barcode_translation.txt.gz"):
    atac, rna = l.strip().split(b"\t")
    atac_to_rna[atac] = rna
assert atac_to_rna[b"atac"] == b"rna"

comp = bytes.maketrans(b'ATGC', b'TACG')

for l in gzip.open(sys.argv[1]):
    fields = l.split(b"\t")
    rev_comp = fields[3].translate(comp)[::-1]
    fields[3] = atac_to_rna[rev_comp]
    sys.stdout.buffer.write(b"\t".join(fields))