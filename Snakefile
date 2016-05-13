# Snakefile

# To run on RCC Midway:
# snakemake -j 100 -c "sbatch --nodes=1 --mem=1000"

import os
from snakemake.utils import R

# Configuration ----------------------------------------------------------------

# Paths must end with forward slash
scratch = "/home/jdblischak/scratch-midway/"
code = "code/"

# Burridge et al 2016
rnaseq_dir = scratch + "burridge2016/"
rnaseq_treatments = ["0um", "1um"]
rnaseq_indivs = ["con1", "con3", "con4", "ch1", "ch3", "ch4"]

for d in [scratch, rnaseq_dir]:
    if not os.path.isdir(d):
        os.mkdir(d)

# Targets ----------------------------------------------------------------------

localrules: process_burridge_2016

rule process_burridge_2016:
    input: expand(rnaseq_dir + "{indiv}-{treatment}.sra", indiv = rnaseq_indivs, treatment = rnaseq_treatments)

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
    input: sra_db = scratch + "SRAmetadb.sqlite"
    output: rnaseq_dir + "{indiv}-{treatment}.sra"
    params: outdir = rnaseq_dir,
            script = code + "download-burridge-2016.R",
            sample = "{indiv}-{treatment}"
    shell: "Rscript {params.script} {input.sra_db} {params.outdir} {params.sample}"
 
