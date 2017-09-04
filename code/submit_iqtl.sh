#!/bin/bash

# Submit per chromosome interaction eQTL scans
# 

if [ -z "$3" ]
then
    echo "Not enough arguments supplied"
    echo "Usage: submit_iqtl.sh normalization permutation_approach cisdist"
    exit 1
fi

for CHROMID in {1..22}; do
  resfile=${DOX_DATA}/sqtl_$1_$2_$3/chr${CHROMID}.txt.gz
  if [ ! -f $resfile ]; then
      #qsub -l nodes=1:ppn=16 -d . -q daglab -v CHROM=chr$CHROMID,NORM=$1,PERM=$2 -N lrt$CHROMID run_iqtl.sh
      sbatch run_iqtl.batch $1 chr$CHROMID $2 $3
  fi
done

resfile=${DOX_DATA}/sqtl_$1_$2_$3/chrX.txt.gz
if [ ! -f $resfile ]; then
    #qsub -l nodes=1:ppn=16 -d . -q daglab -v CHROM=chrX,NORM=$1,PERM=$2 -N lrtX run_iqtl.sh
    sbatch run_iqtl.batch $1 chrX $2 $3
fi
