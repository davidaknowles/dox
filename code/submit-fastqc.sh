#!/bin/bash

# Run from $dox

for FQ in `ls fastq/*`
do
  ID=`basename ${FQ%.fastq.gz}`
  echo "run-fastqc.sh $FQ" | qsub -N fastqc.$ID -S /bin/bash -l walltime=1:00:00 -l nodes=1:ppn=4 -l mem=4gb -e $HOME -V -d . -j eo 
done
