#!/bin/sh

for f in ~/scailscratch/dox/fastq/*.fastq.gz; do
  s=${f%.*}
  outdir=${s%.*}_out/
  echo ${outdir}
  infile=${outdir}Aligned.out.sam
  samtools view -bS $infile | samtools sort - > ${outdir}Aligned_sorted.bam
  rm $infile
done
