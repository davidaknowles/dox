#!/usr/bin/env python3

import glob
import sys
import os
files=glob.glob(os.path.expanduser("~/scailscratch/dox/fastq/*.md5"))
files.sort()
with open("../data/fastq-md5sum-dak.txt","w") as outf:
    for fn in files:
        print(fn)
        with open(fn) as f:
            outf.write( f.readlines()[0].split()[0] + "  " + os.path.splitext( fn.split("/")[-1] )[0] + ".fastq.gz\n" )

