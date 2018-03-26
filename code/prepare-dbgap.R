#!/usr/bin/env Rscript

# Combine the 2 FASTQ files per sample, rename, and record md5 checksum.
#
# Run from project directory

indir <- "/project2/gilad/jdblischak/dox/fastq"
fq <- list.files(indir, "fastq.gz$", full.names = TRUE)
outdir <- "/project2/gilad/jdblischak/dox/fastq-combined"
dir.create(outdir, showWarnings = FALSE)

anno <- read.delim("data/annotation.txt", stringsAsFactors = FALSE)

# Each sample was sequenced on 2 lanes
stopifnot(length(fq) / nrow(anno) == 2)

# Create and submit a shell script for each sample that:
#
# 1. Uncompresses and concatenates the 2 fastq.gz files (use anno$id)
# 2. Compresses the concatenated file (use anno$sample_id)
# 3. Caculates and records the md5 checksum

# Decompress and compress in separate steps to obtain better compression:
#
# -c --stdout --to-stdout
#
# Write output on standard output; keep original files unchanged.  If there are
# several input files, the output consists of a sequence of independently
# compressed members. To obtain better compression, concatenate all input files
# before compressing them.
script <-
"#!/bin/bash
set -eux

zcat %s > %s
gzip %s
md5sum %s > %s
"

for (i in seq_along(anno$id)) {
  if (is.na(anno$sample_id[i])) {
    next
  }
  infiles <- fq[grep(anno$id[i], fq)]
  stopifnot(length(infiles) == 2)
  infiles <- paste(infiles, collapse = " ")
  outfile <- file.path(outdir, paste0(anno$sample_id[i], ".fastq"))
  outfile_gz <- paste0(outfile, ".gz")
  md5file <- file.path(outdir, paste0(anno$sample_id[i], ".md5"))
  script_to_submit <- tempfile(anno$id[i], fileext = ".sbatch")
  writeLines(sprintf(script, infiles, outfile, outfile, outfile_gz,
                     md5file), con = script_to_submit)
  sbatch <- sprintf("sbatch -p %s --mem=%s -o %s %s", "broadwl", "4G",
                    paste0(anno$id[i], ".out"), script_to_submit)
  cat(sbatch, sep = "\n")
  system(sbatch)
  Sys.sleep(1)
}
