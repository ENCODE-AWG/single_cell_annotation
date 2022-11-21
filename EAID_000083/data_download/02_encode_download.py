#!/usr/bin/env python
# Download files from ENCODE related to a given set of experiments
import argparse
import multiprocessing.pool
import os.path
import requests
import json

import tqdm

# Example usage:
# python 03_code/encode_download/02_encode_download.py --experiment_ids 02_config/encode_experiment_ids.txt --cookies 02_config/encode_cookies.json --threads 6 --out_dir 01_raw_data/ENCODE/files

# Description: 
# Downloads all files listed under the given experiment_ids which have a correct output type.
# Saves files under out_dir/ENCFF[...].[extension]
#  (So far, it seems each output_type is guaranteed to always have the same file_type, so I'm 
#   stopping filtering by file_type)

def main():
    parser = argparse.ArgumentParser(description='Download multiome FASTQ files from ENCODE')
    parser.add_argument('--out_dir', type=str, help='base dir for outputs')
    parser.add_argument('--experiment_json', type=str, help='output path for experiment details json')
    parser.add_argument('--output_types', type=str, nargs="+", help="List of output types to download", default=["reads", "index reads"])
    parser.add_argument('--threads', type=int, help="Number of files to download in parallel", default=8)

    args = parser.parse_args()
    
    pool = multiprocessing.pool.ThreadPool(args.threads)    
    os.makedirs(args.out_dir, exist_ok=True)

    experiments = json.load(open(args.experiment_json, "r"))

    files = [f["href"] for e in experiments for f in e["files"] if f["output_type"] in args.output_types]
    file_sizes = {f["href"]: f["file_size"] for e in experiments for f in e["files"] if f["output_type"] in args.output_types}
    total_size = sum(file_sizes.values())
    print(f"Downloading {len(files)} files, {total_size/1e9:.1f}GB total")

    progress = tqdm.tqdm(total=total_size, disable=None, unit_scale=True, unit="B")

    def download_files(file_href):
        out_path = os.path.join(
            args.out_dir,
            os.path.basename(file_href)
        )
        # Check if it's already downloaded and the size is correct
        if os.path.isfile(out_path) and os.stat(out_path).st_size == file_sizes[file_href]:
            progress.update(file_sizes[file_href])
            return
        
        encode_download_file(file_href, out_path, progress)

    
    for _ in pool.imap_unordered(download_files, files):
        pass
        

def encode_download_file(url_path, out_path, progress=None):
    with requests.get(f"https://www.encodeproject.org/{url_path}", stream=True) as r:
        r.raise_for_status()
        with open(out_path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1000000): 
                f.write(chunk)
                if progress:
                    progress.update(len(chunk))


if __name__ == "__main__":
   main()

