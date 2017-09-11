#!/bin/bash

# Submit per chromosome interaction eQTL scans
# 

echo "Data dir: $DOX_DATA"

if [ -z "$3" ]
then
    echo "Not enough arguments supplied"
    echo "Usage: submit_iqtl.sh normalization permutation_approach cisdist"
    exit 1
fi

cisdist=$(printf "%g" "$3")
for CHROMID in {1..22}; do
  resfile=${DOX_DATA}/panama_$1_$2_$cisdist/chr${CHROMID}.txt.gz
  if [ ! -f $resfile ]; then
      #qsub -l nodes=1:ppn=16 -d . -q daglab -v CHROM=chr$CHROMID,NORM=$1,PERM=$2 -N lrt$CHROMID run_iqtl.sh
      echo $resfile not found
      sbatch  -J iqtl_$CHROMID run_iqtl.batch $1 chr$CHROMID $2 $3
  else
      echo $resfile exists
  fi
done

exit 0

resfile=${DOX_DATA}/sqtl_$1_$2/chrX.txt.gz
if [ ! -f $resfile ]; then
    #qsub -l nodes=1:ppn=16 -d . -q daglab -v CHROM=chrX,NORM=$1,PERM=$2 -N lrtX run_iqtl.sh
    sbatch run_iqtl.batch $1 chrX $2 $3
else
    echo $resfile exists
fi
