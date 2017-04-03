#!/bin/bash
set -e

# Submit jobs to convert each chromosome of cgi genotypes to VCF.
#
# Performed on PPS cluster (Sun Grid Engine)

# Directory to read cgi genotype files
dir_cgi=/mnt/gluster/home/jdblischak/ober/imputed-override3
# Directory to write vcf genotype files
dir_vcf=/mnt/gluster/home/jdblischak/ober/vcf
# Path to Plink .fam file
fam=/mnt/gluster/home/jdblischak/ober/Hutterite_imputation/qc.fam
# Path to dox samples file
samples=/mnt/gluster/home/jdblischak/dox/data/samples.txt
# Path to Python conversion script
script=/mnt/gluster/home/jdblischak/dox/code/convert-cgi-to-vcf.py
# Path to save log files
dir_log=$dir_vcf

mkdir -p $dir_vcf $dir_log

for chr in {1..22}
do
  cmd="python $script $fam $samples $dir_cgi/imputed_cgi.chr$chr.tsv.gz $dir_vcf/dox-hg38-chr$chr.vcf.gz"
  echo $cmd | qsub -l h_vmem=8g -V -cwd -N "dox-hg38-chr$chr" -j y -o $dir_log \
                   -l 'hostname=!(bigmem01|bigmem02|spudling90)'
done
