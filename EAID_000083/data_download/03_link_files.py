#!/usr/bin/env python
# Link downloaded fastq files into a directory format usable by cellranger
import argparse
import os.path
import json

import pandas as pd

import tqdm

# symlink files into a more reasonable directory structure (not designed for use
# on the fastqs)
def main():
    parser = argparse.ArgumentParser(description='Download multiome FASTQ files from ENCODE')
    parser.add_argument('--experiment_json', type=str, help="Path of json file with downloaded experiment data")
    parser.add_argument('--experiment_tsv', type=str, help='TSV with columns atac_id, rna_id, and id')
    parser.add_argument('--output_types', type=str, nargs="+", help="List of output types to download", default=["fragments", "unfiltered sparse gene count matrix of unique reads"])
    parser.add_argument('--file_dir', type=str, help="Directory with raw downloads named by ENCODE file accession")
    parser.add_argument('--out_dir', type=str, help="Directory to put output symlinks")

    args = parser.parse_args()
    
    os.makedirs(args.out_dir, exist_ok=True)

    experiments_json = json.load(open(args.experiment_json, "r"))
    
    experiments = pd.read_csv(args.experiment_tsv, sep="\t")
    
    sample_ids = {}
    for e in experiments.itertuples():
        sample_ids[e.atac_id] = e.id
        sample_ids[e.rna_id] = e.id

    files = [(sample_ids[e["accession"]], f) for e in experiments_json for f in e["files"] if f["output_type"] in args.output_types]
    
    for e, f in tqdm.tqdm(files):
        experiment_dir = os.path.join(args.out_dir, e)
        os.makedirs(experiment_dir, exist_ok=True)
        
        src = os.path.join(args.file_dir, os.path.basename(f["href"]))
        dst = os.path.join(experiment_dir, os.path.basename(f["submitted_file_name"]))
        src = os.path.relpath(src, start = os.path.dirname(dst))

        if os.path.islink(dst) and os.readlink(dst) == src:
            pass
        else:
            os.symlink(src, dst)

if __name__ == "__main__":
   main()