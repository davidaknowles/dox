#!/usr/bin/env python3

# Create VCF file of phased genotypes.
#
# Usage:
#
# python convert-cgi-to-vcf.py fam samples cgi vcf
#
# fam - the Plink .fam file with each row corresponding to a column in the cgi genotype file
# samples - a subset of individuals (in 2nd column) from fam to output in VCF
# cgi - gzipped tsv file of genotypes from Complete Genomics
# vcf - gzipped vcf file
#
# Create test file with subset of genotypes:
#
# $ zcat imputed-override3/imputed_cgi.chr21.tsv.gz | head -n 100000 | tail -n 25 | gzip -c > test_cgi.chr21.tsv.gz

# Link to VCF4.3 format specification:
# http://samtools.github.io/hts-specs/VCFv4.3.pdf
# https://en.wikipedia.org/wiki/Variant_Call_Format

# The CGI file encodes each genotype as PXX. X is each allele (0 or 1). P is a
# prefix indicating phasing assignment. The prefix is as followed:
#
# 0 is unphased
# 1 is phased
# 2 is phased with parent of origin
# 3 is imputed using impute2 and phased
# 4 is imputed using impute2 and phased with parent of origin

import datetime
import glob
import gzip
import os
import pandas as pd
from pyliftover import LiftOver
# https://github.com/konstantint/pyliftover
import sys

# Functions --------------------------------------------------------------------

def write_gzip(handle, string):
    # Convert to byte so that string can be written to gzipped file
    #
    # handle - file handle opened with gzip.open
    # string -
    assert isinstance(string, str), "Input must be a string."
    handle.write(string.encode("utf-8"))

def write_vcf_metainfo(handle, version, date, source, reference,
                     contig, phasing):
    write_gzip(handle, "##fileformat=%s\n"%(version))
    write_gzip(handle, "##fileDate=%s\n"%(date))
    write_gzip(handle, "##source=%s\n"%(source))
    write_gzip(handle, "##reference=%s\n"%(reference))
    write_gzip(handle, "##contig=%s\n"%(contig))
    write_gzip(handle, "##phasing=%s\n"%(phasing))
    write_gzip(handle, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n")
    write_gzip(handle, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

def write_vcf_header(handle, individuals):
    write_gzip(handle, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
    for ind in individuals:
        write_gzip(handle, "\t" + ind)
    write_gzip(handle, "\n")

def format_vcf(chr, pos, rsid, phase, allele1, allele2,
               ref, alt, qual, info, format):
    variant_info = "%s\t"*9%(chr, pos, rsid,
                             ref, alt, qual, filter,
                             info, format)
    result = variant_info
    for i in range(len(allele1)):
        if phase[i] == "2" or phase == "4":
            p = "|"
        else:
            p = "/"
        if allele1[i] == "N":
            a1 = "."
        elif allele1[i] in ["0", "1"]:
            a1 = allele1[i]
        else:
            exit("Invalid allele 1 of variant %s: %s"%(rsid,
                                                       allele1[i]))
        if allele2[i] == "N":
            a2 = "."
        elif allele2[i] in ["0", "1"]:
            a2 = allele2[i]
        else:
            exit("Invalid allele 2 of variant %s: %s"%(rsid,
                                                       allele2[i]))

        g = a1 + p + a2
        result = result + g + "\t"

    # Change the final character from a tab to a newline
    result = result[:-1] + "\n"
    return result

# Variables --------------------------------------------------------------------

args = sys.argv[1:]
assert len(args) == 4, \
"""Incorrect input arguments.

Usage:

python convert-cgi-to-vcf.py fam samples cgi vcf

fam - the Plink .fam file with each row corresponding to a column in the cgi genotype file
samples - a subset of individuals (in 2nd column) from fam to output in VCF
cgi - gzipped tsv file of genotypes from Complete Genomics
vcf - gzipped vcf file
"""

# Input CGI file
in_fname = args[2]
#in_fname = "test_cgi.chr21.tsv.gz"
assert os.path.exists(in_fname), "Input cgi file exists."
in_handle = gzip.open(in_fname)

# Output VCF file
out_fname = args[3]
#out_fname = "test_cgi.chr21.vcf.gz"
out_handle = gzip.open(out_fname, "wb")

# Plink .fam file
fam_fname = args[0]
#fam_fname = "qc.fam"
assert os.path.exists(fam_fname), "Input cgiPlink .fam file exists."
fam_handle = open(fam_fname)

# Dox individuals
dox_fname = args[1]
#dox_fname = "../data/samples.txt"
assert os.path.exists(dox_fname), "Input samples file exists."
dox_handle = open(dox_fname)

contig = os.path.basename(in_fname).split(".")[1]

qual = "."
filter = "."
format = "GT"

# liftOver chain
lo = LiftOver("hg19", "hg38")

# Today's date
today = str(datetime.date.today()).replace("-", "")

# Logging number of SNPs
log_variants = 0   # Number of variants
log_snps = 0       # Number of SNPs
log_success = 0    # Number of SNPs successfully converted to hg38
log_nonconvert = 0 # Number of SNPs with no coordinates in hg38
log_multi = 0      # Number of SNPs with multiple coordinates in hg38
log_diff_chr = 0  # Number of SNPs with hg38 coordinates on diff chromosome

# Parse individuals -----------------------------------------------------------

# All individuals in Plink .fam file
ind_all = []
for line in fam_handle:
    cols = line.strip().split()
    ind_all = ind_all + [cols[1]]
assert len(ind_all) == 1415, "There are 1415 individuals in the family file"

# dox individuals
ind_dox = []
for line in dox_handle:
    cols = line.strip().split()
    ind_dox = ind_dox + [cols[1]]
assert len(ind_dox) == 46, "There are 46 individuals in dox project"

assert len(set(ind_all).intersection(ind_dox)) == 46, \
  "dox individuals are a subset of all individuals"

# Sort individuals by their order in the genotype files
ind_dox.sort(key = lambda x: ind_all.index(x))

# Parse file -------------------------------------------------------------------

write_vcf_metainfo(handle = out_handle, version = "VCFv4.3",
                   date = today, source = "CGI",
                   reference = "hg38", contig = contig,
                   phasing = "partial")
write_vcf_header(handle = out_handle, individuals = ind_dox)

for line in in_handle:
    log_variants += 1
    cols = line.decode("utf-8").strip().split("\t")
    id = cols[0]
    chr = cols[1]
    assert chr == contig, "Chromosome of SNP must match contig in filename"
    chr = chr[3:]
    start = cols[2]
    end = cols[3]
    type = cols[4]
    ref = cols[5]
    alt = cols[6]
    if cols[7] == "-":
        rsid = "."
    else:
        # Sometimes it includes multipe dbSNP/rsID combinations, so always
        # take the more recent rsID (farmost right).
        rsid = cols[7].split(":")[-1]
    phase = [x[0] for x in cols[8:]]
    allele1 = [x[1] for x in cols[8:]]
    allele2 = [x[2] for x in cols[8:]]

    # Only report SNPs
    if type != "snp":
        continue
    log_snps += 1

    # Convert hg19 coordinate to hg38 using pyliftover
    # CGI is 0-based, so feed the start position to pyliftover
    lo_result = lo.convert_coordinate("chr" + chr, int(start))

    # Unable to convert from hg19 to hg38
    if len(lo_result) == 0:
        log_nonconvert += 1
        sys.stderr.write("No conversion\t%s\t%s\t%s\n"%(chr, start, rsid))
        continue

    # Multiple coordinates on hg38
    if len(lo_result) > 1:
        log_multi += 1
        sys.stderr.write("Multiple coordinates\t%s\t%s\t%s\n"%(chr, start, rsid))
        continue

    lo_chr, lo_pos = lo_result[0][:2]

    # Different chromosome between hg19 and hg38
    if lo_chr != "chr" + chr:
        log_diff_chr += 1
        sys.stderr.write("Different chromsome\t%s\t%s\t%s\n"%(chr, start, rsid))
        continue

    log_success += 1
    # VCF position is 1-based, so add 1
    pos = str(lo_pos + 1)

    # Subset by individuals in dox project
    phase_dox = [phase[i] for i in range(len(ind_all)) if ind_all[i] in ind_dox]
    allele1_dox = [allele1[i] for i in range(len(ind_all)) if ind_all[i] in ind_dox]
    allele2_dox = [allele2[i] for i in range(len(ind_all)) if ind_all[i] in ind_dox]
    assert len(phase_dox) == len(ind_dox), "One phase prefix per dox ind"
    assert len(allele1_dox) == len(ind_dox), "One allele1 per dox ind"
    assert len(allele2_dox) == len(ind_dox), "One allele2 prefix per dox ind"

    # INFO. Calculate NS - Number of samples with data,
    # i.e. the number of individuals with both alleles known.
    ns = 0
    for aa in zip(allele1_dox, allele2_dox):
        if aa[0] != "N" and aa[1] != "N":
            ns += 1
    assert ns >= 0 and ns <= len(ind_dox), \
      "Reasonable number for number of samples with data"
    info = "NS=%d"%(ns)

    out = format_vcf(chr, pos, rsid, phase_dox,
                     allele1_dox, allele2_dox,
                     ref, alt, qual, info, format)
    write_gzip(out_handle, out)

in_handle.close()
out_handle.close()

# Report logging on SNPs
sys.stderr.write("\n\nLog:\n")
sys.stderr.write("Number of variants:\t%d\n"%(log_variants))
sys.stderr.write("Number of SNPs:\t%d\n"%(log_snps))
sys.stderr.write("Number of SNPs successfully converted to hg38:\t%d\n"%(log_success))
sys.stderr.write("Number of SNPs with no coordinates in hg38:\t%d\n"%(log_nonconvert))
sys.stderr.write("Number of SNPs with multiple coordinates in hg38:\t%d\n"%(log_multi))
sys.stderr.write("Number of SNPs with hg38 coordinates on diff chromosome:\t%d\n"%(log_diff_chr))
