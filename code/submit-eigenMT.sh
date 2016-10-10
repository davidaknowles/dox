#!/bin/bash

for CHROM in {1..22}; do
    outfile=~/scailscratch/dox/eigenMT_chr$CHROM.txt
    if [ ! -f $outfile ]; then
       qsub -d . -q daglab -v CHROM=chr$CHROM -N chr$CHROM run-eigenMT.sh
    fi
done
