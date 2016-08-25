#!/usr/bin/env Rscript

# Download RNA-seq data from Burridge et al., 2016.
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76314
#
# This takes a long time.
#
# Usage:
#
# Rscript download-burridge-2016.R outdir
#
# outdir - Directory to save files, created if it doesn't exist.
#          Warning: Don't use symlinks. sra-tools does not accept.
#
# Steps:
#
# 1. Download SRA database and obtain sample information.
#
# 2. Download and rename SRA files.
#
# 3. Convert to fastq files (using sra-tools).
#
# Results:
#
# 24 gzipped fastq files from paired-end sequencing of 12 samples:
#
# con1-0um
# con1-1um
# con3-0um
# con3-1um
# con4-0um
# con4-1um
# ch1-0um
# ch1-1um
# ch3-0um
# ch3-1um
# ch4-0um
# ch4-1um
#
# Installation:
#
# This R script has the following dependencies:
#
# R packages: dplyr, GEOquery, SRAdb
# sra-tools https://github.com/ncbi/sra-tools/wiki/Downloads
#

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("GEOquery"))
suppressPackageStartupMessages(library("SRAdb"))

# Input command line parameters
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript download-burridge-2016.R outdir")
}
outdir <- args[1]
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

ind <- c("con1", "con3", "con4", "ch1", "ch3", "ch4")
treatment <- c("0um", "1um")
samples <- paste(rep(ind, each = 2), treatment, sep = "-")

# For interactive testing:
# sra_db <- "/scratch/midway/jdblischak/SRAmetadb.sqlite"
# outdir <- "/scratch/midway/jdblischak/"

# 1. Download SRA database and obtain sample information -----------------------

# Download SRAmetadb.sqlite.gz to output directory
sra_db <- file.path(outdir, "SRAmetadb.sqlite")
if (!file.exists(sra_db)) {
  getSRAdbFile(destdir = outdir, method = "wget")
}
stopifnot(file.exists(sra_db))

# Download meta information from GEO
gse_raw <- getGEO("GSE76314")
gse <- as.data.frame(gse_raw[[1]])
sra_info <- gse %>%
  mutate(id = sub(" ", "-", tolower(title)),
         srx = basename(as.character(supplementary_file_1))) %>%
  select(id, srx)
stopifnot(length(unique(sra_info$id)) == length(sra_info$id),
          length(unique(sra_info$srx)) == length(sra_info$srx))

# Confirm that all samples are included
stopifnot(sort(samples) == sort(sra_info$id))

# 2. Download and rename SRA files ---------------------------------------------

# Translate from SRX to SRR
sra_con <- dbConnect(SQLite(), sra_db)
original_name <- character(length = length(samples))
new_name <- character(length = length(samples))
for (srx_index in 1:nrow(sra_info)) {
  srr <- sraConvert(sra_info$srx[srx_index], "run", sra_con)$run
  # They uploaded one run (i.e. fastq) per sample
  stopifnot(length(srr) == 1)
  getSRAfile(srr, sra_con, destDir = outdir, method = "wget")
  original_name[srx_index] <- sprintf("%s/%s.sra", outdir, srr)
  new_name[srx_index] <- sprintf("%s/%s.sra", outdir, sra_info$id[srx_index])
  stopifnot(file.exists(original_name[srx_index]))
  file.rename(original_name[srx_index], new_name[srx_index])
}

dbDisconnect(sra_con)

# 3. Convert to fastq files (using sra-tools) ----------------------------------

for (sra_file in new_name) {
  cmd <- sprintf("fastq-dump --split-files --gzip --outdir %s %s",
                 outdir, sra_file)
  system(cmd)
}

# Confirm final files were created
final_files <- gsub(".sra", "", new_name)
final_files <- paste(final_files, 1:2, sep = "_")
final_files <- paste0(final_files, ".fastq.gz")
stopifnot(file.exists(final_files))
