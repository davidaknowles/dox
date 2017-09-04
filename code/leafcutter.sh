#!/bin/bash                                                                                                                                                                                          
set -e 

BAM_DIR=/afs/cs.stanford.edu/u/davidknowles/scailscratch/dox/bam/

BAM_FILE=dox_bams.txt
ls ${BAM_DIR}*.bam > $BAM_FILE

LEAFCUTTER_DIR=/afs/cs.stanford.edu/u/davidknowles/Dropbox/splicing/leafcutter/

anno=${LEAFCUTTER_DIR}/leafcutter/data/gencode.v26.annotation.gtf.gz

python ${LEAFCUTTER_DIR}example_data/run_sQTL.py -t ${BAM_DIR}../leafcutter_out/ -o dox -d $LEAFCUTTER_DIR -b $BAM_FILE
