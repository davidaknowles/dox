#!/usr/bin/env python3

import os
import shutil

start_dir = "/group/gilad-lab/dox"
dest_dir = "/group/gilad-lab/dox/fastq"

if not os.path.exists(dest_dir):
    os.mkdir(dest_dir)

flow_cell = "C6TJGACXX"

for (path, dirs, files) in os.walk(start_dir):
    if len(files) == 3:
        for f in files:
            if "fastq.gz" in f:
                name_parts = os.path.basename(f).rstrip("fastq.gz").split("_")
                id = "%02d"%(int(name_parts[0]))
                lane = name_parts[2][0] + name_parts[2][3]
                new_name = "-".join([id, flow_cell, lane, name_parts[1]])
                new_name = new_name + ".fastq.gz"
                print(new_name)
                shutil.copyfile(path + "/" + f, dest_dir + "/" + new_name)

