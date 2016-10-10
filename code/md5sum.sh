#!/bin/sh
#$ -S /bin/bash
cd $PBS_O_WORKDIR

md5sum $FILE > $OUT
