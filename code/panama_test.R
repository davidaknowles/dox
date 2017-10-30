
library("dplyr")
library("tidyr")
library(data.table)
source("utils.R")

require(rstan)
source("load_data.R")

require(suez)

geneloc=read.table(paste0(DATADIR,"genelocGRCh38.txt.gz"),header=T,stringsAsFactors = F)

if (interactive()) {
  chrom="chr15"
  normalization_approach="qq"
  permutation_approach="boot"
  cisdist=1e4
  # registerDoMC(16)
} else {
  ca=commandArgs(trailingOnly = T)
  chrom=ca[2]
  permutation_approach=ca[3]
  normalization_approach=ca[1]
  cisdist=as.numeric(ca[4])
  registerDoMC(16)
}

geneloc=geneloc[geneloc$chr==chrom,]
snploc=snploc[snploc$chr==chrom,]

if (interactive()) {
   geneloc=head(geneloc,2)
}

stopifnot(all(as.character(snploc$snpid) %in% rownames(genotype) ))
genotype=genotype[as.character(snploc$snpid),]

sample_kernel=readRDS("../data/Kern.rds")

resdir=paste0(DATADIR,"panama_",normalization_approach,"_",permutation_approach,"_",cisdist,"/")
dir.create(resdir)

checkpoint_dir=paste0(resdir,chrom,"_checkpoint/")
dir.create(checkpoint_dir)

anno_for_suez=anno %>% select(individual=dbgap, condition=conc)

results = map_interaction_qtl(input, genotype, geneloc, snploc, anno_for_suez, sample_kernel, normalization_approach, permutation_approach, cisdist, checkpoint_dir, debug=F) 

if (!interactive()){
   gz1 = gzfile(paste0(resdir,chrom,".txt.gz"),"w")
   results %>% format(digits=5) %>% write.table(gz1, quote = F, row.names = F, col.names = T, sep="\t")
   close(gz1)
}
