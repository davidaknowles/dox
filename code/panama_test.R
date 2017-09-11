
library("dplyr")
library("tidyr")
library(data.table)
source("utils.R")

require(rstan)
source("map_interaction_qtl.R")
source("load_data.R")

geneloc=read.table(paste0(DATADIR,"genelocGRCh38.txt.gz"),header=T,stringsAsFactors = F)

if (interactive()) {
  chrom="chr15"
  normalization_approach="qq"
  permuted="boot"
  cisdist=1e5
  # registerDoMC(16)
} else {
  ca=commandArgs(trailingOnly = T)
  chrom=ca[2]
  permuted=ca[3]
  normalization_approach=ca[1]
  cisdist=as.numeric(ca[4])
  registerDoMC(16)
}

geneloc=geneloc[geneloc$chr==chrom,]
snploc=snploc[snploc$chr==chrom,]

stopifnot(all(as.character(snploc$snpid) %in% rownames(genotype) ))
genotype=genotype[as.character(snploc$snpid),]

#input=remove_PCs(input, num_PCs_to_remove)

sample_kernel=readRDS("../data/Kern.rds")

resdir=paste0(DATADIR,"panama_",normalization_approach,"_",permuted,"_",cisdist,"/")
dir.create(resdir)

checkpoint_dir=paste0(resdir,chrom,"_checkpoint/")
dir.create(checkpoint_dir)

results = map_interaction_qtl(input, geneloc, snploc, anno, sample_kernel, normalization_approach, checkpoint_dir)

gz1 = gzfile(paste0(resdir,chrom,".txt.gz"),"w")
results %>% format(digits=5) %>% write.table(gz1, quote = F, row.names = F, col.names = T, sep="\t")
close(gz1)
