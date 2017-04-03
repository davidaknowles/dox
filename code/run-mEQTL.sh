#!/bin/sh
#$ -S /bin/bash

hostname

echo NPC: $NPC
echo CONC: $CONC

set -e

Rscript run_MatrixEQTL.R $CONC $NPC

# Rscript collate-mEQTL-qsub.R $CONC $NPC
