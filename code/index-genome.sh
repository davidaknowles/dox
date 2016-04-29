#!/bin/bash
set -e

# Build index for mapping with Subread

# First argument is output directory
OUT_DIR=$1

mkdir -p $OUT_DIR
cd $OUT_DIR

echo "Build index"
subread-buildindex -o hg19 *.fa
