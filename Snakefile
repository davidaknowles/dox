# Snakefile

# To run on RCC Midway:
# snakemake -j 100 --cluster-config config-rcc.json -c "sbatch --mem={cluster.mem} --nodes={cluster.n} --tasks-per-node={cluster.tasks}"

import os
from snakemake.utils import R

# Configuration ----------------------------------------------------------------

# Paths must end with forward slash
scratch = "/scratch/midway/jdblischak/"

code = "code/"
data = "data/"

# Burridge et al 2016
rnaseq_dir = scratch + "burridge2016/"
rnaseq_treatments = ["0um", "1um"]
rnaseq_indivs = ["con1", "con3", "con4", "ch1", "ch3", "ch4"]

for d in [scratch, rnaseq_dir]:
    if not os.path.isdir(d):
        os.mkdir(d)

# Targets ----------------------------------------------------------------------

localrules: process_burridge_2016, prepare_kallisto

rule process_burridge_2016:
    input: expand(rnaseq_dir + "{indiv}-{treatment}/abundance.tsv", \
                  indiv = rnaseq_indivs, treatment = rnaseq_treatments)

rule prepare_kallisto:
    input: scratch + "transcriptome-ensembl-GRCh38.idx"

# Rules ------------------------------------------------------------------------

rule download_sra_meta:
    output: temp(scratch + "SRAmetadb.sqlite")
    params: destdir = scratch
    run:
        R("""
        library("SRAdb")
        getSRAdbFile(destdir = "{params.destdir}", method = "wget")
        """)

rule download_burridge_rnaseq:
    input: sra_db = scratch + "SRAmetadb.sqlite"
    output: temp(rnaseq_dir + "{indiv}-{treatment}.sra")
    params: outdir = rnaseq_dir,
            script = code + "download-burridge-2016.R",
            sample = "{indiv}-{treatment}"
    shell: "Rscript {params.script} {input.sra_db} {params.outdir} {params.sample}"
 
rule convert_sra_to_fastq:
    input: "{sample}.sra"
    output: "{sample}_1.fastq.gz", "{sample}_2.fastq.gz"
    shell: "fastq-dump --split-files --gzip --outdir `dirname {input}` {input}"

rule download_transcriptome:
    output: scratch + "transcriptome-ensembl-GRCh38.fa.gz"
    shell: "wget -O {output} http://bio.math.berkeley.edu/kallisto/transcriptomes/Homo_sapiens.GRCh38.rel79.cdna.all.fa.gz"

rule kallisto_index:
    input:  scratch + "transcriptome-ensembl-GRCh38.fa.gz"
    output: scratch + "transcriptome-ensembl-GRCh38.idx"
    shell: "kallisto index -i {output} {input}"

rule kallisto_quant:
    input: read1 = "{sample}_1.fastq.gz", read2 = "{sample}_2.fastq.gz",
           index = scratch + "transcriptome-ensembl-GRCh38.idx"
    output: "{sample}/abundance.tsv"
    params: outdir = "{sample}", threads = 8, bootstraps = 100
    shell: "kallisto quant -i {input.index} -o {params.outdir} -t {params.threads} -b {params.bootstraps} {input.read1} {input.read2}"

rule kallisto_collate:
    input: expand(rnaseq_dir + "{indiv}-{treatment}/abundance.tsv", \
                  indiv = rnaseq_indivs, treatment = rnaseq_treatments)
    output: eff_counts = data + "burridge-2016-eff-counts.txt",
            tpm =  data + "burridge-2016-tpm.txt",
    params: script = code + "kallisto-collate.R"
    shell: "Rscript {params.script} {output.eff_counts} {output.tpm} {input}"
