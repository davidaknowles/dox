# Snakefile

# http://jdblischak.github.io/singleCellSeq/analysis/islam2014.html

import os
from snakemake.utils import R

# Configuration ----------------------------------------------------------------

# Paths must end with forward slash
scratch = "/home/jdblischak/scratch-midway/"

# Burridge et al 2016
rnaseq_treatments = ["0um", "1um"]
rnaseq_samples = ["con1", "con3", "con4", "ch1", "ch3", "ch4"]

# Targets ----------------------------------------------------------------------



# Rules ------------------------------------------------------------------------

rule download_sra_meta:
    output: scratch + "SRAmetadb.sqlite"
    params: destdir = scratch
    run:
        R("""
        library("SRAdb")
        getSRAdbFile(destdir = "{params.destdir}", method = "wget")
        """)

rule download_burridge_rnaseq:
    input: scratch + "SRAmetadb.sqlite"
    output: expand(scratch + "{sample}-{treatment}.fastq.gz", sample = rnaseq_samples, treatment = rnaseq_treatments)
    params: outdir = scratch
    run:
        R("""
        library("dplyr")
        library("GEOquery")
        library("SRAdb")

        # Download meta information from GEO
        gse <- getGEO("GSE76314")
        gse <- as.data.frame(gse[[1]], stringsAsFactors = FALSE)
        sra_info <- gse %>%
                    mutate(id = sub(" ", "-", tolower(title)),
                           srx = substr(supplementary_file_1,
                                        start = 77, stop = 85)) %>%
                    select(id, srx)

        # Translate from SRX to SRR
        sra_con <- dbConnect(SQLite(), "{input}")
        for (srx_index in 1:nrow(sra_info)) {{
          srr <- sraConvert(sra_info$srx[srx_index], "run", sra_con)$run
          for (srr_index in 1:length(srr)) {{
            original_name <- sprintf("%s/%s_%d.fastq.gz", "{params.outdir}", srr[srr_index], srr_index)
            new_name <- sprintf("%s/%s-%d.fastq.gz", "{params.outdir}",
                                sra_info$id[srx_index], srr_index)
            getSRAfile(srr[srr_index], sra_con, destDir = "{params.outdir}",
                       fileType = "fastq", method = "wget")
            stopifnot(file.exists(original_name))
            file.rename(original_name, new_name)
          }}
        }}
        """)
