#!/usr/bin/env Rscript

# Download RNA-seq data from Burridge et al., 2016.
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76314
#
# Usage:
#
# Rscript download-burridge-2016.R sra_db outdir [samples...]
#
# sra_db - SRAmetadb.sqlite file
# outdir - directory to save SRA files, created if it doesn't exist
# samples - (optional) specific samples to download. Default is all 12 samples.
#
# Options for [samples...]:
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

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("GEOquery"))
suppressPackageStartupMessages(library("SRAdb"))

# Input command line parameters
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Must supply at least two arguments.")
}
sra_db <- args[1]
outdir <- args[2]

if (length(args) > 2) {
  samples <- args[3:length(args)]
} else {
  ind <- c("con1", "con3", "con4", "ch1", "ch3", "ch4")
  treatment <- c("0um", "1um")
  samples <- paste(rep(ind, each = 2), treatment, sep = "-")
}

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
         srx = basename(as.character(supplementary_file_1))) %>%
  select(id, srx)
stopifnot(length(unique(sra_info$id)) == length(sra_info$id),
          length(unique(sra_info$srx)) == length(sra_info$srx))

# Filter based on input samples
sra_info <- sra_info %>% filter(id %in% samples)
if (nrow(sra_info) < 1) {
  stop(sprintf("None of the supplied samples were valid:\n%s", samples))
}

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
