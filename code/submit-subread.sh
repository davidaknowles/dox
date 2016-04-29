#!/bin/bash

# Run from $dox

for FQ in `ls fastq/*`
do
  ID=`basename ${FQ%.fastq.gz}`
  echo "run-subread.sh $FQ" | qsub -N subread.$ID -S /bin/bash -l walltime=1:00:00 -l nodes=1:ppn=4 -l mem=12gb -o $HOME -e $HOME -V -d .
done
