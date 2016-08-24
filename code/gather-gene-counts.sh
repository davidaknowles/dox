#!/bin/bash

~/Dropbox/dox/code/gather-gene-counts.py ~/scailscratch/dox/counts/s*.genecounts.txt | gzip > ~/Dropbox/dox/data/gene-counts-round-two.txt.gz
