# Note: This script can in principle be run end-to-end, but
# is better used by copy-pasting commands one at a time with any modifications
# needed for number of threads/cores to use

mkdir -p data
mkdir -p plots
mkdir -p data/azimuth_references/azimuth_panc

# Download Azimuth reference datasets
wget https://zenodo.org/record/4546926/files/idx.annoy?download=1 -O data/azimuth_references/azimuth_panc/idx.annoy
wget https://zenodo.org/record/4546926/files/ref.Rds?download=1 -O data/azimuth_references/azimuth_panc/ref.Rds


# Download gene annotation gtf
curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz > \
    data/gencode.v29.annotation.gtf.gz

# Download Tabula Sapiens Pancreas reference
wget https://figshare.com/ndownloader/files/34702003 -O data/TS_Pancreas.h5ad.zip
unzip data/TS_Pancreas.h5ad.zip -d data


# Fetch metadata for desired experiment IDs
python data_download/01_download_metadata_json.py \
    --experiment_tsv config/panc_samples.tsv \
    --out_json config/panc_experiments.json \
    --threads 25

# Download experiments (Just fragments and gene count matrices)
python data_download/02_encode_download.py \
    --experiment_json config/panc_experiments.json \
    --out_dir data/ENCODE/files \
    --threads 25 \
    --output_types "fragments" "unfiltered sparse gene count matrix of unique reads"

# Construct a nice directory listing by sample
python data_download/03_link_files.py \
    --experiment_json config/panc_experiments.json \
    --experiment_tsv config/panc_samples.tsv \
    --file_dir data/ENCODE/files \
    --out_dir data/ENCODE/samples \
    --output_types "fragments" "unfiltered sparse gene count matrix of unique reads"

# Unpack the tar.gz files into directly readable files
snakemake -j 16 -s data_download/04_unpack_tar.smk --config experiments_tsv="config/panc_samples.tsv"

# Run data import tasks with 16 cores
snakemake -j 16 -s data_import/04_data_import.smk

# Run all the analysis scripts in order
Rscript analysis/00a_export_qc_data.R
Rscript analysis/00b_raw_rna_and_qc_plots.R
Rscript analysis/01a_integrate_tabula_sapiens.R
Rscript analysis/01b_integrate_samples.R
Rscript analysis/02_export_cell_coords.R
Rscript analysis/03_export_cell_type_calls_and_plots.R
Rscript analysis/04_export_misc.R