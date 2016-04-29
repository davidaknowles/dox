#!/bin/bash

# Run from $dox

for BAM in `ls bam/*bam`
do
  ID=`basename ${BAM%.bam}`
  echo "run-samtools.sh $BAM" | qsub -N samtools.$ID -S /bin/bash -l walltime=1:00:00 -l nodes=1:ppn=4 -l mem=8gb -o $HOME -e $HOME -V -d .
done
