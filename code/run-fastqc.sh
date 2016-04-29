#!/bin/bash
set -e

FILE=$1
BASE=`basename ${FILE%.fastq.gz}`
OUTDIR=fastqc

mkdir -p $OUTDIR

if [ ! -s $FILE ]
then
  echo "File is missing or empty: $FILE"
  exit 65
fi

if [ -s $OUTDIR/$BASE.raw.count.txt ]
then
  echo "Output file already exists: $OUTDIR/$BASE.raw.count.txt"
  exit 64
fi

# Unzip file (fastqc throws error when passed unzipped file via process substitution)
zcat $FILE > $OUTDIR/$BASE.fastq

# Run FastQC
fastqc $OUTDIR/$BASE.fastq

# Count number of reads
echo $(( $(cat $OUTDIR/$BASE.fastq | wc -l) / 4 )) > $OUTDIR/$BASE.count.txt

# Remove unzipped fastq file
rm $OUTDIR/$BASE.fastq
