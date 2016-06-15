#!/bin/bash

# To be run from directory with genotype files.

original=qc
new=dox
samples=samples.txt

# Filter individuals
# https://www.cog-genomics.org/plink2/filter#indiv
plink --bfile $original --make-bed --keep $samples --out $new

# Convert to VCF
# https://www.cog-genomics.org/plink2/data#recode
plink --bfile $new --recode vcf-iid --out $new
