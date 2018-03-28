#!/bin/sh

start=`date +%s`
zcat < $1 | awk '(NR==1) || (($3 >= -100000) && ($3 <= 100000))' | gzip > $2
end=`date +%s`

runtime=$((end-start))
echo Runtime: $runtime

