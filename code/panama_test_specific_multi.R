

library("dplyr")
library("tidyr")
library(data.table)
source("utils.R")
#registerDoMC(7)
require(rstan)

source("load_data.R")

normalization_approach="qq"
permuted="boot"

cisdist=1e6

errorhandling='remove'

K=readRDS("../data/Kern.rds")

mega_res=foreach(rsid=gwas_ld_filtered$other_snp, .combine = bind_rows) %dopar% {
  gwas_hit=snploc %>% filter(RSID==rsid)
  if (nrow(gwas_hit) != 1) return(NULL)

sample_kernel=readRDS("../data/Kern.rds")

resdir=paste0(DATADIR,"panama_",normalization_approach,"_",permuted,"_",cisdist,"/")
dir.create(resdir)

checkpoint_dir=paste0(resdir,chrom,"_checkpoint/")
dir.create(checkpoint_dir)

results = map_interaction_qtl(input, geneloc, snploc, anno, sample_kernel, normalization_approach, cisdist, checkpoint_dir)

gz1 = gzfile(paste0(resdir,chrom,".txt.gz"),"w")
results %>% format(digits=5) %>% write.table(gz1, quote = F, row.names = F, col.names = T, sep="\t")
close(gz1)