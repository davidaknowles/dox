#!/usr/bin/env python3

# Purpose: Confirm that individuals with genotyped parents have more
# SNPs that are assigned a parent of origin compared to individuals
# with 0 or 1 genotyped parents.
#
# Each genotype for each individual is specified with 3
# characters. The first is a number from 0 to 4 which specieis the
# results of the phasing. The interpretation of the prefix is
# below. The next two numbers specify the allele of each chromosome. 0
# is the reference allele, 1 the alternative, and N is unknown.
#
# The prefix is as followed,
# 0 is unphased.
# 1 is phased
# 2 is phased with parent of origin
# 3 is imputed using impute2 and phased
# 4 is imputed using impute2 and phased with parent of origin
#
# It computes the following statistics for each individual across all
# the automsomes:
#
# total - number of SNPs
# unknown - number of unknown SNPs (any SNP w/ at least one N is
#           counted here and not any of the categories below)
# p0 - number of unphased SNPs
# p1 - number of phased SNPs
# p2 - number of SNPs with parent of origin assignments
# p3 - number of phased SNPs which were imputed with impute2
# p4 - number of SNPs which were imputed with impute2
#
# Note that these categories are mutually exclusive, and thus the sum
# of unknown, p0, p1, p2, p3, and p4 will equal total for each
# individual.
#
# Usage:
#
# python3 count-phasing.py
#
# Warning: This takes hours to complete.
#
# Output:
#
# Creates the tab-separated file count-phasing.txt. It has 7 columns
# which correspond to the statistic calculated above. The order of the
# rows (individuals) is the same as in the PLINK fam file.
#
# $ head count-phasing.txt 
# p0	p1	p2	p3	p4	total	unknown
# 104823	2056671	4160619	958129	367073	9188660	1541345
# 24133	0	8507556	189622	68104	9188660	399245
# 42052	0	8201111	361407	141569	9188660	442521
# 35970	0	8302581	336859	86860	9188660	426390
# 25547	0	8463194	201701	93787	9188660	404431
# 51044	969213	6980968	407545	160286	9188660	619604
# 31479	772247	7649489	284528	80345	9188660	370572
# 11050	1151108	7763124	0	0	9188660	263378
# 56707	0	7606230	531308	158003	9188660	836412

import glob
import gzip
import os
import pandas as pd
import sys

# Path to genotype files
target_dir = "imputed-override3/"

chromosomes = [str(i) for i in range(1, 22 + 1)]
# For testing:
#chromosomes = ["21"]

geno_total = [0] * 1415
geno_unknown = [0] * 1415
geno_0 = [0] * 1415
geno_1 = [0] * 1415
geno_2 = [0] * 1415
geno_3 = [0] * 1415
geno_4 = [0] * 1415

for chr in chromosomes:
    sys.stderr.write("\nProcessing chromsome %s\n\n"%(chr))
    chr_fname = target_dir + "imputed_cgi.chr%s.tsv.gz"%(chr)
    sys.stderr.write("Filename: %s\n"%(chr_fname))
    chr_handle = gzip.open(chr_fname, "r")
    for line in chr_handle:
        cols = line.decode("utf-8").strip().split()
#        import pdb; pdb.set_trace()
        if cols[4] != "snp":
            continue
        for i in range(1415):
            geno = cols[8 + i]
            geno_total[i] += 1
            if geno[1:] == "NN":
                geno_unknown[i] += 1
            elif geno[0] == "0":
                geno_0[i] += 1
            elif geno[0] == "1":
                geno_1[i] += 1
            elif geno[0] == "2":
                geno_2[i] += 1
            elif geno[0] == "3":
                geno_3[i] += 1
            elif geno[0] == "4":
                geno_4[i] += 1
            elif geno[0] == "5":
                geno_5[i] += 1
            else:
                sys.stderr.write("Unknonwn genotype classification for %s"%(geno))
    chr_handle.close()

# Output results
out_fname = "count-phasing.txt"
df = pd.DataFrame({
    "total": geno_total,
    "unknown": geno_unknown,
    "p0": geno_0,
    "p1": geno_1,
    "p2": geno_2,
    "p3": geno_3,
    "p4": geno_4
})
df.to_csv(out_fname, sep = "\t", index = False)
