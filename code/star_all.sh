#!/bin/bash

THREADS=16

for f in ~/scailscratch/dox/fastq/*.fastq.gz; do
  s=${f%.*}
  outdir=${s%.*}_out/
  if [ ! -f ${outdir}Log.final.out ]; then
      mkdir $outdir
      bn=$( basename $outdir )
      echo $bn
      rm -r ${outdir}*
      qsub -l nodes=1:ppn=$THREADS -d . -q daglab -v OUTPREFIX=$outdir,FILE=$f,THREADS=$THREADS -N $bn star.sh
  fi
  #infile=${outdir}Aligned.out.sam
  #qsub -l nodes=1:ppn=$THREADS -d . -q daglab -v outdir=$outdir,infile=$infile -N $bn star.sh
done
