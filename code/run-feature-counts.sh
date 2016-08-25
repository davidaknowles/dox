#!/bin/sh
#$ -S /bin/bash

hostname

#cd $PBS_O_WORKDIR
mkdir /tmp/featureCount/
cd /tmp/featureCount/

echo Input: $FILE
echo Output: $OUTFILE

set -e

EXONS=~/Dropbox/dox/data/exons_GRCh38.saf

echo "Counting reads per gene..."
featureCounts -T 8 -a $EXONS -F SAF -o $OUTFILE $FILE
