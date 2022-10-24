
# Command-line arguments:
# 1. sample name
# 2. path of input arrow file
# 3. path of output base directory

# example: Rscript create_project.R sample1.arrow my/directory/sample1
# Outputs to my/directory/sample1/ArchR_Project_unfiltered 

#args <- c("W73_PANC", "04_data/cellranger_arc_v2/workspace/ENCSR229VVY/outs/atac_fragments.tsv.gz", "04_data/archr/W73_PANC")
args <- commandArgs(trailing=TRUE)
stopifnot(length(args) == 3)
sample_name <- args[1]
input_fragments <- args[2]
output_base <- args[3]

# wd <- getwd()
# setwd("/oak/stanford/groups/wjg/bparks/multiome/")
# source("renv/activate.R")
# setwd(wd)

suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
})

dir.create(output_base, recursive=TRUE, showWarnings=FALSE)

addArchRThreads(1)

base_dir <- getwd()

arrow_output <- file.path(base_dir, output_base, paste0(sample_name, ".arrow"))
qc_output <- file.path(base_dir, output_base, "QualityControl")
logs_output <- file.path(base_dir, output_base, "ArchRLogs")


tmp <- tempdir()

# Copy arrow files into tmp for potentially better I/O performance
tmp_fragments <- file.path(tmp, "fragments.tsv.gz")
tmp_arrow_output <- file.path(tmp, sample_name)
tmp_qc_output <- file.path(tmp, "QualityControl")

file.copy(file.path(base_dir, input_fragments), tmp_fragments)
file.copy(file.path(base_dir, paste0(input_fragments, ".tbi")), paste0(tmp_fragments, ".tbi"))

setwd(tmp)

addArchRGenome("hg38")

ArrowFiles <- createArrowFiles(
  inputFiles =  tmp_fragments,
  sampleNames = sample_name,
  outputNames = file.path(tmp, sample_name),
  QCDir = file.path(tmp, "QualityControl"),
  minTSS=0,
  minFrags=0,
  addTileMat = FALSE,
  addGeneScoreMat = FALSE
)

# Copy outputs into the destination locations
file.copy(ArrowFiles, arrow_output)

dir.create(qc_output)
dir.create(logs_output)
file.copy(
  file.path(tmp, "QualityControl", sample_name), 
  qc_output,
  recursive=TRUE
)
file.copy(
  file.path(tmp, "ArchRLogs"),
  logs_output,
  recursive=TRUE
)
