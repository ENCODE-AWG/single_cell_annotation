import collections
import json
import os.path
import re

import pandas as pd

samples = pd.read_csv("config/panc_samples.tsv", sep="\t")

sample_ids = list(samples.id)


rule all:
    input: 
        expand("data/archr/{sample_id}/ArchR_Project_unfiltered", sample_id=sample_ids),
        expand("data/seurat/{sample_id}/rna_raw.rds", sample_id=sample_ids),
        "data/TS_Panc.h5seurat"
        
rule rna_mats:
    input: "data/ENCODE/samples/{sample_id}/rna"
    output: "data/seurat/{sample_id}/rna_raw.rds"
    params:
        script = srcdir("02_rna_mat_to_rds.R")
    log: "data/seurat/logs/import_{sample_id}.log"
    shell: "Rscript {params.script} {input} {output} > {log} 2> {log}"

rule archr_import:
    input: "data/ENCODE/samples/{sample_id}/fragments.tsv.gz"
    output: directory("data/archr/{sample_id}/ArchR_Project_unfiltered")
    params:
        output_base = "data/archr/{sample_id}/ArchR_Project_unfiltered",
        script = srcdir("01_archr_create_project.R")
    log: "data/archr/logs/import_{sample_id}.log"
    shell: "Rscript {params.script} {wildcards.sample_id} {input} {params.output_base} > {log} 2> {log}"

rule convert_tabula_sapiens:
    input: "data/TS_Panc.h5ad"
    output: "data/TS_Panc.h5seurat"
    params:
        script = srcdir("03_convert_h5seurat.R")
    shell: "Rscript {params.script} {input} {output}"