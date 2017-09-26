#!/bin/bash

THREADS=16

for f in ~/scailscratch/dox/Burridge/*_1.fastq.gz; do
  s=${f%.*}
  s=${s%.*}
  s=${s::${#s}-2}
  echo $s
  outdir=${s}_out/
  if [ ! -f ${outdir}Log.final.out ]; then
      mkdir $outdir
      bn=$( basename $outdir )
      echo $bn
      f1=${s}_1.fastq.gz
      f2=${s}_2.fastq.gz
      rm -r ${outdir}*
      qsub -l nodes=1:ppn=$THREADS -d . -q daglab -v OUTPREFIX=$outdir,FILE1=$f1,FILE2=$f2,THREADS=$THREADS -N $bn star_paired.sh
  fi
  #infile=${outdir}Aligned.out.sam
  #qsub -l nodes=1:ppn=$THREADS -d . -q daglab -v outdir=$outdir,infile=$infile -N $bn star.sh
done
