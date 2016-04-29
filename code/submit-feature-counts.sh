#!/bin/bash

# Run from $dox

for BAM in `ls bam-processed/*bam`
do
  ID=`basename ${BAM%.bam}`
  echo "run-feature-counts.sh $BAM" | qsub -N feature.counts.$ID -S /bin/bash -l walltime=1:00:00 -l nodes=1:ppn=4 -l mem=4gb -e $HOME -V -d . -j eo
done
