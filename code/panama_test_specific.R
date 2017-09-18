
library("dplyr")
library("tidyr")
library(data.table)
source("utils.R")

require(rstan)
source("map_interaction_qtl.R")
source("load_data.R")

geneloc=read.table(paste0(DATADIR,"genelocGRCh38.txt.gz"),header=T,stringsAsFactors = F)

normalization_approach="qq"
permutation_approach="boot"
cisdist=1e6

stopifnot(all(as.character(snploc$snpid) %in% rownames(genotype) ))

sample_kernel=readRDS("../data/Kern.rds")

resdir=paste0(DATADIR,"rs11855704_",normalization_approach,"_",permutation_approach,"_",cisdist,"/")
dir.create(resdir)

checkpoint_dir=paste0(resdir,chrom,"_checkpoint/")
dir.create(checkpoint_dir)

snploc=snploc %>% filter(RSID=="rs11855704")

results = map_interaction_qtl(input, genotype, geneloc, snploc, anno, sample_kernel, normalization_approach, permutation_approach, cisdist, checkpoint_dir) 

