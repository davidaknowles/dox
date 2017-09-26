#!/bin/bash

for f in ~/scailscratch/dox/Burridge/*_1.fastq.gz; do
  s=${f%.*}
  s=${s%.*}
  stub=${s::${#s}-2}
  outbam=${stub}_out/Aligned.sortedByCoord.out.bam
  if [ -f $outbam ]; then
      bn=$( basename $stub )
      mv $outbam ~/scailscratch/dox/Burridge/bam/${bn}.bam
  fi
done
