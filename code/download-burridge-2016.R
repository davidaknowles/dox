#!/usr/bin/env Rscript

# Download RNA-seq data from Burridge et al., 2016.
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76314
#
# Usage
#
# Rscript download-burridge-2016.R sra_db outdir
#
#


suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("GEOquery"))
suppressPackageStartupMessages(library("SRAdb"))

# Input command line parameters
args <- commandArgs(trailingOnly = TRUE)
sra_db <- args[1]
outdir <- args[2]

# For interactive testing
# sra_db <- "/home/jdblischak/scratch-midway/SRAmetadb.sqlite"
# outdir <- "/home/jdblischak/scratch-midway/"

stopifnot(file.exists(sra_db))
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# Download meta information from GEO
gse_raw <- getGEO("GSE76314")
gse <- as.data.frame(gse_raw[[1]])
sra_info <- gse %>%
  mutate(id = sub(" ", "-", tolower(title)),
         srx = substr(supplementary_file_1,
                      start = 77, stop = 86)) %>%
  select(id, srx)
stopifnot(length(unique(sra_info$id)) == length(sra_info$id),
          length(unique(sra_info$srx)) == length(sra_info$srx))

# Translate from SRX to SRR
sra_con <- dbConnect(SQLite(), sra_db)
for (srx_index in 1:nrow(sra_info)) {
  srr <- sraConvert(sra_info$srx[srx_index], "run", sra_con)$run
  # They uploaded one run (i.e. fastq) per sample
  stopifnot(length(srr) == 1)
  getSRAfile(srr, sra_con, destDir = outdir, method = "wget")
  original_name <- sprintf("%s/%s.sra", outdir, srr)
  new_name <- sprintf("%s/%s.sra", outdir, sra_info$id[srx_index])
  stopifnot(file.exists(original_name))
  file.rename(original_name, new_name)
}

dbDisconnect(sra_con)
