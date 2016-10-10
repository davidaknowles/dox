#!/bin/sh
#$ -S /bin/bash

ggcd $PBS_O_WORKDIR

hostname

echo $CHROM

Rscript run_iqtl.R $CHROM