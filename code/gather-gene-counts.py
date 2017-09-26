#!/usr/bin/env python3

# Gathers all the counts per gene.
# usage: gather-gene-counts.py [files] > gene-counts.txt
#   e.g. gather-gene-counts.py counts/*genecounts.txt > gene-counts.txt
# Should be run from data directory.

import glob
import sys
import gzip 

if len(sys.argv) > 1:
    files = sys.argv[1:]
else:
    files = glob.glob("counts/*genecounts.txt")

# print(len(files))

# Create header
sys.stdout.write("filename\tindividual\tflow_cell\tlane\tindex\tconc")
# Get gene names from first file
gene_list = []
f = files[0]
handle = gzip.open(f,"r") if f[-2:]=="gz" else open(f, "r")
for line in handle:
    line=str(line)
    if line[0] == "#" or line[:6] == "Geneid":
        continue
    #print(line, file=sys.stderr)
    cols = line.strip("\n").split("\t")
    gene = cols[0]
    gene_list.append(gene)
    sys.stdout.write("\t" + gene)
handle.close()
sys.stdout.write("\n")

# Process each file
for f in files:
    # Get meta data from filename
    fname = f.split("/")[-1]
    stub=fname.strip(".genecounts.txt")
    fname_parts = stub.split("-")
    #print fname_parts
    individual, flow_cell, conc, index, lane = fname_parts
    # Update variable names
    individual = "i" + individual

    sys.stdout.write(stub + "\t" + individual + "\t" + flow_cell + "\t" + \
                     lane + "\t" + index + "\t" + conc)

    # Get counts from f
    g = 0 # iterator for indexing gene names
    handle = gzip.open(f,"r") if f[-2:]=="gz" else open(f, "r")
    for line in handle:
        if line[0] == "#" or line[:6] == "Geneid":
            continue
        cols = line.strip("\n").split("\t")
        gene = cols[0]
        assert gene == gene_list[g], "Gene names are not in correct order in file: %s"%(f)
        g += 1
        counts = cols[6]
        sys.stdout.write("\t" + counts)
    handle.close()
    sys.stdout.write("\n")
