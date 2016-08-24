#!/usr/bin/env python

def get_fn(fn):
    with open(fn) as f:
        a=[ g.strip().split()[1] for g in f.readlines() ]
    return a
    
a=get_fn("../data/fastq-md5sum.txt")
b=get_fn("../data/fastq-md5sum-dak.txt") 

missing=list(set(a)-set(b))

missing.sort()

fn="../data/missing-fastq.txt"

if len(missing)==0:
    print("No missing files")
else:
    print("%i missing files, writing to %s" % (len(missing), fn))
    with open(fn,"w") as f:
        f.writelines(missing)

import os
print("Problematic md5sums:")
os.system("diff ../data/fastq-md5sum.txt ../data/fastq-md5sum-dak.txt")
