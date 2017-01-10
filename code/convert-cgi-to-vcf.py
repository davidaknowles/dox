#!/usr/bin/env python3

# Create VCF file of phased genotypes.

# To do:
#
# * Filter by call rate
# * Filter individuals
# * Convert to hg38
#
# Create test file:
#
# $ zcat imputed-override3/imputed_cgi.chr21.tsv.gz | head -n 100000 | tail -n 25 | cut -f1-20 | gzip -c > test_cgi.chr21.tsv.gz

import glob
import gzip
import os
import pandas as pd

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
    write_gzip(handle, "##INFO=<ID=CR,Number=1,Type=Float,Description=\"Call Rate\">\n")
    write_gzip(handle, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

def write_vcf_header(handle, individuals):
    write_gzip(handle, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
    for ind in individuals:
        write_gzip(handle, "\t" + ind)
    write_gzip(handle, "\n")

def format_vcf(chr, start, rsid, phase,
               ref, alt, qual, info, format):
    variant_info = "%s\t"*9%(chr, start, rsid,
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

in_fname = "test_cgi.chr21.tsv.gz"
in_handle = gzip.open(in_fname)
out_fname = "test_cgi.chr21.vcf.gz"
out_handle = gzip.open(out_fname, "wb")

contig = os.path.basename(in_fname).split(".")[1]

qual = "."
filter = "."

write_vcf_metainfo(handle = out_handle, version = "VCFv4.2",
                   date = "20170110", source = "CGI",
                   reference = "hg19", contig = contig,
                   phasing = "partial")
#write_vcf_header(handle = outfile, individuals = x.individual)

for line in in_handle:
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
        dbsnp = "."
        rsid = "."
    else:
        # To do: sometimes it includes multipe dbSNP/rsID combinations
        database = cols[7].split(":")
        dbsnp = database[0]
        rsid = database[1]
    phase = [x[0] for x in cols[8:]]
    allele1 = [x[1] for x in cols[8:]]
    allele2 = [x[2] for x in cols[8:]]
    individual = ["ind" + str(i) for i in range(len(allele1))]

    # Need to add call rate to INFO
    info = ""
    format = "GT"

    out = format_vcf(chr, start, rsid, phase,
                     ref, alt, qual, info, format)
    write_gzip(out_handle, out)
        
in_handle.close()
out_handle.close()
