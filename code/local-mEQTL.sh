#!/bin/bash

for concid in {1..5}; do
    for pcs in {0,1,2,5,10,15,20}; do
      N=${concid}_$pcs
      echo $N
      Rscript run_MatrixEQTL.R $concid $pcs
    done
done

