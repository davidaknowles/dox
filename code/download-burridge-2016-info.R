#!/usr/bin/env Rscript

# Download SRA information for RNA-seq data from Burridge et al., 2016.
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76314
#
# The output file is burridge-sra-info.txt.
#
# Warning: This downloads a 20G file that can be deleted afterwards.
#
# Usage:
#
# Rscript download-burridge-2016-info.R outdir
#
# outdir - Directory to save files, created if it doesn't exist.
#
# Steps:
#
# 1. Download SRA database.
#
# 2. Download sample information from GEO.
#
# 3. Obtain URLs to download SRA files from FTP site.
#
# Results:
#
# Writes tab-delimited file burridge-sra-info.txt to outdir.
#
# Also saves the 20G SRA sqlite database to outdir.
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
  stop("Usage: Rscript download-burridge-2016-info.R outdir")
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
outdir <- "/scratch/midway/jdblischak/"

# 1. Download SRA database and obtain sample information -----------------------

# Download SRAmetadb.sqlite.gz to output directory
sra_db <- file.path(outdir, "SRAmetadb.sqlite")
if (!file.exists(sra_db)) {
  getSRAdbFile(destdir = outdir, method = "wget")
}
stopifnot(file.exists(sra_db))

# 2. Download sample information from GEO ----------------------------------------

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

# 3. Obtain URLs to download SRA files from FTP site ---------------------------

# Translate from SRX to SRR
sra_con <- dbConnect(SQLite(), sra_db)
original_name <- character(length = length(samples))
new_name <- character(length = length(samples))
ftp <- character(length = length(samples))
for (srx_index in 1:nrow(sra_info)) {
  srr <- sraConvert(sra_info$srx[srx_index], "run", sra_con)$run
  # They uploaded one run (i.e. fastq) per sample
  stopifnot(length(srr) == 1)
  sra_ftp_info <- listSRAfile(srr, sra_con)
  ftp[srx_index] <- sra_ftp_info$ftp
  original_name[srx_index] <- paste0(srr, ".sra")
  new_name[srx_index] <- paste0(sra_info$id[srx_index], ".sra")
}

dbDisconnect(sra_con)

# Save SRA information

write.table(data.frame(original_name, new_name, ftp),
            file = file.path(outdir, "burridge-sra-info.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)
