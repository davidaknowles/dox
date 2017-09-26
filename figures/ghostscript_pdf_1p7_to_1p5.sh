#!/bin/bash 
set -e
outfile=$2
if [ "$#" -eq 1 ]; then
    outfile=$1
fi
gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dCompatibilityLevel=1.5 -sOutputFile=temp.pdf $1
# so you can run on the same file
mv temp.pdf $outfile