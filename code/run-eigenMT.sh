#!/bin/sh
#$ -S /bin/bash

ggcd $PBS_O_WORKDIR

hostname

echo $CHROM

~/myvirtualenv/bin/python ../eigenMT/eigenMT.py --CHROM $CHROM --QTL ~/scailscratch/dox/mEQTL_results/results0_PC10.txt.gz --GEN ~/scailscratch/dox/genotype.txt --GENPOS ~/scailscratch/dox/snploc.txt --PHEPOS ~/scailscratch/dox/genelocGRCh38.txt --OUT ~/scailscratch/dox/eigenMT_$CHROM.txt
