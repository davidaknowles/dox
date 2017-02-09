#!/bin/bash
set -e

# Submit jobs to check genotypes to reference genome.
#
# Performed on PPS cluster (Sun Grid Engine)

# Directory to read vcf genotype files
dir_vcf=/mnt/gluster/home/jdblischak/ober/vcf
# Path to R script
script=/mnt/gluster/home/jdblischak/dox/code/convert-cgi-to-vcf-test.R
# Path to save log files
dir_log=$dir_vcf

mkdir -p $dir_vcf $dir_log

for chr in {1..22}
do
  cmd="Rscript $script $dir_vcf/dox-hg38-chr$chr.vcf.gz"
  echo $cmd | qsub -l h_vmem=12g -V -cwd -N "dox-hg38-chr$chr-test" -j y -o $dir_log \
                   -l 'hostname=!(bigmem01|bigmem02|spudling90)'
done
