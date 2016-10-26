#!/bin/bash

# Submit per chromosome interaction eQTL scans

for CHROMID in {1..22}; do
  resfile=~/scailscratch/dox/lrt_imp_chr${CHROMID}.RData
  if [ ! -f $resfile ]; then
      qsub -l nodes=1:ppn=16 -d . -q daglab -v CHROM=chr$CHROMID -N lrt$CHROMID run_iqtl.sh
  fi
done

qsub -l nodes=1:ppn=16 -d . -q daglab -v CHROM=chrX -N lrtX run_iqtl.sh
