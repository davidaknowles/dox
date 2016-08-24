#!/bin/sh
#$ -S /bin/bash
cd $PBS_O_WORKDIR

echo File: $FILE
echo Out: $OUTPREFIX
echo Threads: $THREADS

STAR --runThreadN $THREADS --genomeDir ~/scailscratch/STAR_index/ --readFilesIn $FILE --readFilesCommand zcat --outFileNamePrefix $OUTPREFIX --outSAMtype BAM SortedByCoordinate

#samtools view -bS $infile | samtools sort - > ${outdir}Aligned_sorted.bam
#rm $infile
