#!/usr/bin/env python
# Download JSON metadata corresponding to the ENCODE experiments
import argparse
import multiprocessing.pool

import requests
import json

import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Download multiome FASTQ files from ENCODE')
    parser.add_argument('--experiment_tsv', type=str, help='TSV with columns atac_id and rna_id')
    parser.add_argument('--out_json', type=str, help='Output path for experiment json objects')
    parser.add_argument('--threads', type=int, help="Number of files to download in parallel", default=8)

    args = parser.parse_args()
    
    experiments = pd.read_csv(args.experiment_tsv, sep="\t")
    experiment_ids = [x for x in list(experiments.atac_id) + list(experiments.rna_id) if x.startswith("ENCSR")]

    pool = multiprocessing.pool.ThreadPool(args.threads)    
    def download_experiment_metadata(experiment_id):
        return encode_fetch_json(f"/experiments/{experiment_id}")
    
    experiments = pool.map(download_experiment_metadata, experiment_ids)
    json.dump(experiments, open(args.out_json, "w"), indent=2, sort_keys=True)

def encode_fetch_json(path):
    headers = {'accept': 'application/json'}
    response = requests.get(f"https://www.encodeproject.org/{path}", headers=headers)
    return response.json()

if __name__ == "__main__":
   main()