#!/bin/bash

# Submit per chromosome interaction eQTL scans

for CHROMID in {1..22}; do
  resfile=~/dagscratch/dox/panama_$1_$2/chr${CHROMID}.txt.gz
  if [ ! -f $resfile ]; then
      qsub -l nodes=1:ppn=16 -d . -q daglab -v CHROM=chr$CHROMID,NORM=$1,PERM=$2 -N lrt$CHROMID run_iqtl.sh
  fi
done

resfile=~/dagscratch/dox/panama_$1_$2/chrX.txt.gz
if [ ! -f $resfile ]; then
    qsub -l nodes=1:ppn=16 -d . -q daglab -v CHROM=chrX,NORM=$1,PERM=$2 -N lrtX run_iqtl.sh
fi
