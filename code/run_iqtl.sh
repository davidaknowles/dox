#!/bin/sh
#$ -S /bin/bash

cd $PBS_O_WORKDIR

echo Chromosome: $CHROM
echo Normalization: $NORM
echo Perm: $PERM

# Rscript run_iqtl_lrt.R $NORM $CHROM
Rscript panama_test.R $NORM $CHROM $PERM
