#!/bin/sh
#$ -S /bin/bash
cd $PBS_O_WORKDIR

echo File: $FILE1 $FILE2
echo Out: $OUTPREFIX
echo Threads: $THREADS

hostname

STAR --runThreadN $THREADS --genomeDir ~/scailscratch/STAR_index/ --readFilesIn $FILE1 $FILE2 --readFilesCommand zcat --outFileNamePrefix $OUTPREFIX --outSAMtype BAM SortedByCoordinate

#samtools view -bS $infile | samtools sort - > ${outdir}Aligned_sorted.bam
#rm $infile
