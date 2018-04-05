#!/bin/bash

# Submit per chromosome interaction eQTL scans
# 

echo "Data dir: $DOX_DATA"

for CHROMID in {1..22}; do
    sbatch run_iqtl.batch chr$CHROMID
done

exit 0
