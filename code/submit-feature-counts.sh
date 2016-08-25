#!/bin/bash

OUTDIR=~/scailscratch/dox/counts
mkdir -p $OUTDIR

for f in ~/scailscratch/dox/bam/*.bam; do
  s=${f%.*}
  bn=$( basename $s )
  echo $bn $f

  BASE=`basename ${f%.bam}`
  
  outfile=$OUTDIR/$BASE.genecounts.txt

  if [ ! -f $outfile ]; then
      qsub -l nodes=1:ppn=16 -d . -q daglab -v FILE=$f,OUTFILE=$outfile -N $bn run-feature-counts.sh
  fi
done
