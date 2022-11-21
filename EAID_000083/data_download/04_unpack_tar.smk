import json

import pandas as pd

# Must have experiments_tsv set in config
samples = pd.read_csv(config["experiments_tsv"], sep="\t")

rule all:
    input: 
        expand("data/ENCODE/samples/{sample}/rna", sample=samples.id),
        expand("data/ENCODE/samples/{sample}/fragments.tsv.gz", sample=samples.id)

translate_py = srcdir("04b_translate_barcodes.py")

def translate_barcodes_command(w, output):
    sample = [e for e in samples.itertuples() if e.id == w.sample][0]
    dir = f"data/ENCODE/samples/{w.sample}"
    if not str(sample.multiome_series).startswith("ENCSR"):
        # No translation needed, so just copy files
        return (
            f"cp {dir}/encode_scatac_dcc_2/results/*/fragments/fragments.tsv.gz {output.fragments}; "
            f"cp {dir}/encode_scatac_dcc_2/results/*/fragments/fragments.tsv.gz.tbi {output.idx}"
        )
    else:
        # Translate barcodes and re-compress
        return (
            f"python {translate_py} {dir}/encode_scatac_dcc_2/results/*/fragments/fragments.tsv.gz | "
            f"bgzip > {output.fragments}; "
            f"tabix --preset bed {output.fragments}"
        )

rule unpack_atac:
    input: 
        tar = "data/ENCODE/samples/{sample}/fragments.tar.gz",
        barcode_translation = "config/10x_barcode_translation.txt.gz"
    output: 
        fragments = "data/ENCODE/samples/{sample}/fragments.tsv.gz",
        idx = "data/ENCODE/samples/{sample}/fragments.tsv.gz.tbi"
    params:
        dir = "data/ENCODE/samples/{sample}",
        translate_command = translate_barcodes_command
    shell: "tar -xzf {input.tar} --directory {params.dir}; "
           "{params.translate_command}; "
           "rm -r {params.dir}/encode_scatac_dcc_2/; "

rule unpack_rna:
    input:
        tar = "data/ENCODE/samples/{sample}/GeneFull_Ex50pAS_Unique_raw.tar.gz"
    output:
        matrix = directory("data/ENCODE/samples/{sample}/rna")
    params:
        dir = "data/ENCODE/samples/{sample}"
    shell: "tar -xzf {input.tar} --directory {params.dir}; "
           "mv {params.dir}/GeneFull_Ex50pAS/raw {params.dir}/rna; "
           "gzip {params.dir}/rna/*; "
           "rm -r {params.dir}/GeneFull_Ex50pAS/; "
