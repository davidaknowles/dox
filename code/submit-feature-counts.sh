#!/bin/bash

# OUTDIR=~/scailscratch/dox/counts
# INDIR=~/scailscratch/dox/bam/
INDIR=$1
OUTDIR=$2

mkdir -p $OUTDIR

for f in ${INDIR}*.bam; do
  s=${f%.*}
  bn=$( basename $s )
  echo $bn $f

  BASE=`basename ${f%.bam}`
  
  outfile=$OUTDIR/$BASE.genecounts.txt

  if [ ! -f $outfile ]; then
      qsub -l nodes=1:ppn=16 -d . -q daglab -v FILE=$f,OUTFILE=$outfile -N $bn run-feature-counts.sh
  fi
done
