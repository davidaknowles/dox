#!/bin/bash
set -e

# Download hg19.

# First argument is output directory
OUT_DIR=$1

mkdir -p $OUT_DIR
cd $OUT_DIR

echo "Download hg19 autosomes"
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr{1..22}.fa.gz

echo "Download hg19 X, Y, and M"
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr{X,Y,M}.fa.gz

echo "Unzip fasta files"
gunzip chr*fa.gz
