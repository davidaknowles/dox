#!/bin/bash

OUTDIR=~/scailscratch/dox/counts
EXONS=~/Dropbox/dox/data/exons_GRCh38.saf

mkdir -p $OUTDIR

for f in ~/scailscratch/dox/bam/*.bam; do
  s=${f%.*}
  bn=$( basename $s )
  echo $bn $f

  BASE=`basename ${f%.bam}`
  
  outfile=$OUTDIR/$BASE.genecounts.txt

  if [ ! -f $outfile ]; then
      qsub -l nodes=1:ppn=1 -d . -q daglab -v FILE=$f,OUTFILE=$outfile -N $bn run-feature-counts-dak.sh
  fi
done
#!/bin/sh
#$ -S /bin/bash

FILE=$1
OUTFILE=$2

echo Input: $FILE
echo Output: $OUTFILE

set -e



echo "Counting reads per gene..."
featureCounts -a $EXONS -F SAF -R -o $OUTFILE $FILE
featureCounts -T 48 -a ~/Dropbox/dox/data/exons_GRCh38.saf -F SAF -o ~/scailscratch/dox/counts/all.txt ~/scailscratch/dox/bam/*.bam
g
