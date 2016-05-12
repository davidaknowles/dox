# Snakefile

# http://jdblischak.github.io/singleCellSeq/analysis/islam2014.html

import os
from snakemake.utils import R

# Configuration ----------------------------------------------------------------

# Paths must end with forward slash
scratch = "/home/jdblischak/scratch-midway/"
code = "code/"

# Burridge et al 2016
rnaseq_dir = scratch + "burridge2016/"
rnaseq_treatments = ["0um", "1um"]
rnaseq_samples = ["con1", "con3", "con4", "ch1", "ch3", "ch4"]

for d in [scratch, rnaseq_dir]:
    if not os.path.isdir(d):
        os.mkdir(d)

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
    input: sra_db = scratch + "SRAmetadb.sqlite",
           script = code + "download-burridge-2016.R"
    output: expand(rnaseq_dir + "{sample}-{treatment}.sra", sample = rnaseq_samples, treatment = rnaseq_treatments)
    params: outdir = rnaseq_dir
    shell: "Rscript {input.script} {input.sra_db} {params.outdir}"
 
