#!/bin/bash

for f in ~/scailscratch/dox/fastq/*.fastq.gz; do
  s=${f%.*}
  stub=${s%.*}
  outbam=${stub}_out/Aligned.sortedByCoord.out.bam
  if [ -f $outbam ]; then
      bn=$( basename $stub )
      mv $outbam ~/scailscratch/dox/bam/${bn}.bam
  fi
done
