#!/usr/bin/env Rscript

# Download RNA-seq data from Burridge et al., 2016.
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76314
#
# This takes a long time.
#
# Usage:
#
# Rscript download-burridge-2016.R outdir ftpfile
#
# outdir - Directory to save files, created if it doesn't exist.
#          Warning: Don't use symlinks. sra-tools does not accept.
#
# ftpfile - tab-delimited file with FTP URL for download (column ftp) and names
#           of files (column new_name)
#
# Steps:
#
# 1. Download SRA files
#
# 2. Convert to fastq files (using sra-tools)
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

# Input command line parameters
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript download-burridge-2016.R outdir ftpfile")
}
outdir <- args[1]
ftpfile <- args[2]

# For interactive testing:
outdir <- "/scratch/midway/jdblischak/"
ftpfile <- "../data/burridge-sra-info.txt"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

ftp <- read.delim(ftpfile, stringsAsFactors = FALSE)
stopifnot(nrow(ftp) > 0,
          c("new_name", "ftp") %in% colnames(ftp))

# 1. Download SRA files --------------------------------------------------------

for (i in 1:nrow(ftp)){
  download.file(url = ftp$ftp[i],
                destfile = file.path(outdir, ftp$new_name[i]),
                method = "wget")
}

# 2. Convert to fastq files (using sra-tools) ----------------------------------

for (i in 1:nrow(ftp)){
  cmd <- sprintf("fastq-dump --split-files --gzip --outdir %s %s",
                 outdir, file.path(outdir, ftp$new_name[i]))
  system(cmd)
}

# Confirm final files were created
final_files <- gsub(".sra", "", file.path(outdir, ftp$new_name))
final_files <- paste(final_files, 1:2, sep = "_")
final_files <- paste0(final_files, ".fastq.gz")
stopifnot(file.exists(final_files))
