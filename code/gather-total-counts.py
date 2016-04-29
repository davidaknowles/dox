#!/usr/bin/env python3

"""
Gathers all the counts contained in the summary files.

usage: gather-total-counts.py > total-counts.txt

Should be run from data directory.
"""

import glob
import sys

files = sorted(glob.glob("fastqc/*count.txt")) + \
        sorted(glob.glob("bam-processed/*count.txt")) + \
        sorted(glob.glob("counts/*summary"))

# print(len(files))

# Output header:
#   stage - the stage of the processing pipeline
sys.stdout.write("stage\tindividual\tflow_cell\tlane\tindex\tcounts\n")

# Function definitions ---------------------------------------------------------

def read_lines(f):
    """
    Input: Path to file
    Output: Lines of the file (list)
    """
    handle = open(f, "r")
    lines = handle.readlines()
    handle.close()
    return lines

def read_count_file(f):
    """
    Input: Path to file
    Output: The number contained in the file (str)
    Explanation: The count file only contains one number, the number of
                 sequences at that stage in the processing pipeline.
    """
    lines = read_lines(f)
    assert "Count file has only one line", len(lines) == 1
    counts = lines[0].strip("\n")
    return counts

def read_featureCounts_summary(f):
    """
    Input: Path to file
    Output: The number of sequences labeled Assigned (str)
    Explanation: featureCounts outputs a file with the extension .summmary that
                 details the number of sequences per result category. The
                 category Assigned is for sequences that map unambiguously to a
                 gene.
    """
    assert "featureCounts summary file has correct extension", \
           f[-8:] == ".summary"
    lines = read_lines(f)
    assert "featureCounts summary file has 12 lines", \
           len(lines) == 12
    assert "The Assigned category is the first entry after the header", \
           lines[1].split("\t")[0] == "Assigned"
    counts = lines[1].strip("\n").split("\t")[1]
    return counts

# Process each file ------------------------------------------------------------

for f in files:

    dir, fname = f.split("/")

    if dir == "fastqc":
        stage = "raw"
        counts = read_count_file(f)
    elif dir == "bam-processed":
        stage = "mapped to genome"
        counts = read_count_file(f)
    elif dir == "counts":
        stage = "mapped to exons"
        counts = read_featureCounts_summary(f)

    # Get meta data from filename
    fname_parts = fname.split(".")[0].split("-")
    individual, flow_cell, lane, index = fname_parts

    # Update variable names
    individual = "i" + individual

    sys.stdout.write(stage + "\t" + individual + "\t" + flow_cell + "\t" + \
                     lane + "\t" + index + "\t" + counts + "\n")
