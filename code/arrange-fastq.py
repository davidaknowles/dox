#!/usr/bin/env python3

# Move and rename fastq files downloaded from FGF's FTP site.
#
# Usage:
#
# python3 arrange-fastq.py indir outdir
#
# indir - Highly nested directory structure from downloading data
# outdir - New output directory (created if does not exist)
#
# Ex:
#
# python3 arrange-fastq.py fgfftp.uchicago.edu/Genomics_Data/NGS/160520_K00242_0070_BHCMNYBBXX-YG-Dox-781112/FastQ fastq
#
# Explanation:
#
# The core provides the samples with the following naming scheme:
#
# 160520_K00242_0070_BHCMNYBBXX-YG-Dox-781112/FastQ/YG-Dox-p11-s105-c29-1-5000_S25_L003_R1_001.fastq.gz
#
# To be extracted are the following variables:
#
# sample number: s105
# cell line num: c29-1
# treatment concentration: 5000
# flow cell id: HCMNYBBXX (the leading A or B is discarded)
# lane: L003
#
# These are converted into the following file naming scheme:
#
# s105-c29.1-5.000-HCMNYBBXX-l3.fastq.gz
#
# sample number: s105 (always has three digits)
# cell line num: c29.1
# treatment concentration: 5.000
# flow cell id: HCMNYBBXX
# lane: l3
#

import glob
import os
import shutil
import sys

# Input arguments
args = sys.argv
assert len(args) == 3, "Incorrect number of arguments.\nUsage: python3 arrange-fastq.py indir outdir"
indir = args[1]
outdir = args[2]
assert os.path.exists(indir), "Input directory does not exist: %s"%(indir)
if not os.path.exists(outdir):
    os.mkdir(outdir)
# Add final forward slash if necessary
if indir[-1] != "/":
    indir = indir + "/"
if outdir[-1] != "/":
    outdir = outdir + "/"

# Obtain file names
files = glob.glob(indir + "/*fastq.gz")[:]

# Rename and move files
for f in files:
    path = f.rstrip('fastq.gz').split('/')
    flow_cell = path[-3].split("_")[-1].split("-")[0][1:]
    file_parts = path[-1].split('_')[:-1]
    lane = "l" + file_parts[2][-1]
    if file_parts[0] == "Undetermined":
        sample_name = file_parts[0].lower()
    else:
        name_parts = file_parts[0].split("-")
        sample_num = name_parts[3]
        sample_num = "s%03d"%(int(sample_num[1:]))
        if len(name_parts) == 6:
            cell_num = name_parts[4]
        elif len(name_parts) == 7:
            cell_num = name_parts[4] + "." + name_parts[5]
        else:
            sys.exit("Input file naming scheme has changed. Code must be updated.")
        treatment = name_parts[-1]
        treatment = treatment[0] + "." + treatment[1:]
        sample_name = "-".join([sample_num, cell_num, treatment])
    new_name = outdir + sample_name + '-' + flow_cell + "-" + lane + '.fastq.gz'
    sys.stderr.write("Moving:\n%s\n%s\n\n"%(new_name, f))
    shutil.move(f, new_name)
