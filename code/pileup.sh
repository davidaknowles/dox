#!/bin/sh
#$ -S /bin/bash

# Print the parameters
echo In 1: $INFILEA
echo In 2: $INFILEB
echo Out: $OUTFILE

cd $PBS_O_WORKDIR

# Stop on error
set -e

# I ended up using conda to install samtools on the cluster :/ 
export PATH=/afs/cs.stanford.edu/u/davidknowles/scailscratch/miniconda2/bin:$PATH

# mpileup: -B prevents error calculation, 
# -d1000000 sets the maximum depth to record
# -l specifies the regions (SNPs) to quantify
# pileup_filter 2 2 is a little prog I wrote to process the mpileup output (should prob use vcf/bcf output option from mpileup instead)
samtools mpileup -B -d1000000 -l ~/scailscratch/Homo_sapiens/NCBI/GRCh38/snp147.txt -f ~/scailscratch/Homo_sapiens/NCBI/GRCh38/Sequence/Bowtie2Index/genome.fa $INFILEA $INFILEB  | ~/Dropbox/eagle/ptsd/pileup_filter 2 2 | gzip > $OUTFILE
