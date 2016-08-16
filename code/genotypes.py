#!/usr/bin/env python3

import glob
import os
import pandas as pd
import subprocess as sp
from pyliftover import LiftOver

# Set working directory
work_dir = "/mnt/gluster/home/jdblischak/ober/"

# The information on the Hutterite genotypes are spread across
# multiple different files. The files are in various formats, and even
# more complicated, each file contains different sets of individuals
# and genotypes.


# File descriptions and paths ----------------------------------------

# The binary plink files contain the genotypes for all 1,415
# individuals for 7,989,987 variants. The only filtering performed was
# on the call rate.

plink_original_prefix = work_dir + "Hutterite_imputation/qc"
# The .bed file is binary and contains the genotype calls.
plink_bed_fname = plink_original_prefix + ".bed"
# The .fam file is plain text, and contains one line per individual
plink_fam_fname =  plink_original_prefix + ".fam"
# The .bim file is plain text, and contains one line per variant
plink_bim_fname =  plink_original_prefix + ".bim"

# Because the plink files were only filtered based on a simple call
# rate filter, they contain variants that are not polymorphic in this
# population, and also indels. The variant information file contains
# more details on each variant. Furthermore, it has already removed
# the non-polymorphic variants.
variant_info_fname = work_dir + "Hutterite_imputation/qc.imputed_cgi.annovar_plink_annotations.hg19_multianno.txt.regulomeDB.txt"

# The phasing information is contained in a file format from Complete
# Genomics. It contains the information on all 12,341,920 SNPs. There
# is one file per chromosome, and each column is an indivudual. The
# ordering of the columns in these files are in the same order as the
# accompanying tfam file.
phasing_glob_pattern = ""
phasing_fam_fname = ""

# The plink individual IDs are the findiv IDs included in the
# annotation file.
anno_fname = "/mnt/gluster/home/jdblischak/dox/data/annotation.txt"

# Filter individuals and genotypes -----------------------------------

# Individuals are filtered with the plink flag --keep. It accepts a
# tab-separated file with the family ID in the first column and the
# inidividaul ID in the second column. There is no header row. The
# family ID is HUTTERITES. Below the findiv IDs are pulled from the
# annotation file to create the file for filtering individuals.
# https://www.cog-genomics.org/plink2/filter#indiv
keep_ind_fname = work_dir + "Hutterite_imputation/samples.txt"
keep_ind_handle = open(keep_ind_fname, "w")
anno_handle = open(anno_fname, "r")
anno_header = anno_handle.readline()
anno_header_cols = anno_header.strip().split("\t")
assert anno_header_cols[3] == "findiv", \
    "The findiv ID needs to be in the 4th column of %s"%(anno_fname)
last_findiv = ""
for line in anno_handle:
    cols = line.strip().split("\t")
    findiv = cols[3]
    # The annotation file has one line per treatment, so the findiv is
    # present a total of 5 times
    if findiv == last_findiv:
        continue
    keep_ind_handle.write("HUTTERITES\t%s\n"%(findiv))
    last_findiv = findiv

anno_handle.close()
keep_ind_handle.close()

# Variants are filtered with the plink flag --extract. It accepts a
# file with one variant ID per line. The variants to be included will
# be those that are designed as "snp" in the 47th column of the
# variant information file.
# https://www.cog-genomics.org/plink2/filter#snp
extract_var_fname = work_dir + "Hutterite_imputation/variants.txt"
extract_var_handle = open(extract_var_fname, "w")
variant_info_handle = open(variant_info_fname, "r")
variant_info_header = variant_info_handle.readline()
variant_info_cols = variant_info_header.strip().split("\t")
assert variant_info_cols[45] == "ImputationID", \
    "The 46th column need to be ImputationID in file %s"%(variant_info_fname)
assert variant_info_cols[46] == "Variant.type", \
    "The 47th column need to be Variant.type in file %s"%(variant_info_fname)
for line in variant_info_handle:
    cols = line.strip().split("\t")
    variant_id = cols[45]
    variant_type = cols[46]
    if variant_type == "snp":
        extract_var_handle.write("%s\n"%(variant_id))

variant_info_handle.close()
extract_var_handle.close()

# Filter individuals and variants using plink
plink_filtered_prefix = work_dir + "Hutterite_imputation/dox"
plink_cmd_filter = "plink --bfile %s --make-bed --keep %s --extract %s --out %s"%(\
    plink_original_prefix, keep_ind_fname, extract_var_fname,
    plink_filtered_prefix)
sp.call(plink_cmd_filter, shell = True)

# Convert to VCF -----------------------------------------------------

# Convert to VCF format using --recode.
# https://www.cog-genomics.org/plink2/data#recode
# https://samtools.github.io/hts-specs/VCFv4.2.pdf
vcf_fname = work_dir + "Hutterite_imputation/dox.vcf"
plink_cmd_vcf = "plink --bfile %s --recode vcf-iid --out %s"%(\
    plink_filtered_prefix, plink_filtered_prefix)
sp.call(plink_cmd_vcf, shell = True)

# Update coordinates to hg38 -----------------------------------------

# Using pyliftover to avoid having to convert to temporary BED file.
# https://github.com/konstantint/pyliftover

lo = LiftOver('hg19', 'hg38')
vcf_hg38_fname = work_dir + "Hutterite_imputation/dox-hg38.vcf"

vcf_hg19_handle = open(vcf_fname, "r")
vcf_hg38_handle = open(vcf_hg38_fname, "w")

for line in vcf_hg19_handle:
    if line[0] == "#":
        vcf_hg38_handle.write(line)
        continue
    cols = line.strip().split("\t")
    # VCF is 1-based, BED is 0-based with the first inclusive and
    # second exclusive. Thus need to subtract one from the start
    # coordinate.
    chrom_hg19 = "chr" + cols[0]
    pos_hg19 = int(cols[1]) - 1
    lo = LiftOver("hg19", "hg38")
    result = lo.convert_coordinate(chrom_hg19, pos_hg19)
    chrom_hg38 = result[0][0]
    pos_hg38 = result[0][1] + 1
    assert chrom_hg19 == chrom_hg38, \
           "hg19 and hg38 chromosomes should match."
    line_hg38 = chrom_hg38 + "\t" + str(pos_hg38) + "\t" + "\t".join(cols[2:]) + "\n"
    vcf_hg38_handle.write(line_hg38)

vcf_hg19_handle.close()
vcf_hg38_handle.close()
