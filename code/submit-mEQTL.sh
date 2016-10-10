#!/bin/bash
concs=(0 0.625 1.25 2.5 5)
for concid in {0..4}; do
    conc=${concs[$concid]}
    # echo $conc
    for pcs in {0,1,2,5,10,15,20}; do
      N=${conc}_$pcs
      echo $N
      ofile=~/scailscratch/dox/mEQTL_results/summary${conc}_PC${pcs}.txt.gz
      if [ ! -f $ofile ]; then
        echo $ofile
        qsub -d . -q daglab -v CONC=$conc,NPC=$pcs -N $N run-mEQTL.sh
      fi
    done
done

