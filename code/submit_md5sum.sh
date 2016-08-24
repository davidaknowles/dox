#!/bin/bash

THREADS=16

for f in ~/scailscratch/dox/fastq/*.fastq.gz; do
  s=${f%.*}
  out=${s%.*}.md5
  if [ ! -f $out ]; then
      bn=$( basename $out )
      echo $bn
      qsub -d . -q daglab -v OUT=$out,FILE=$f -N $bn md5sum.sh
  fi
done
