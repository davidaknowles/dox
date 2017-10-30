

library("dplyr")
library("tidyr")
require(magrittr)
library(data.table)
require(stringr)
source("utils.R")

require(suez)

#  "chr5:102961229:102990272:clu_17500"
if (interactive()) {
  chrom="chr5"
  normalization_approach="none"
  permuted="boot"
  cisdist=1e6
} else {
  registerDoMC(16)
  ca=commandArgs(trailingOnly = T)
  chrom=ca[2]
  permuted=ca[3]
  normalization_approach=ca[1]
  cisdist=1e5
}

source("load_data.R")

snploc=snploc[snploc$chr==chrom,]

stopifnot(all(as.character(snploc$snpid) %in% rownames(genotype) ))
genotype=genotype[as.character(snploc$snpid),]

input <- read.table(paste0(DATADIR,"leafcutter_qqnorm.txt.gz"), header=T, sep="\t", check.names = F)
anno=str_split_fixed(colnames(input), "_", 2) %>% 
  as.data.frame(stringsAsFactors=F) %>%
  set_colnames(c("findiv","conc"))

geneloc=str_split_fixed(rownames(input),":",4) %>% 
  as.data.frame(stringsAsFactors=F) %>%
  set_colnames(c("chr","left","right","clu")) %>%
  mutate(left=as.numeric(left), right=as.numeric(right),geneid=rownames(input))

geneloc=geneloc %>% filter(chr==chrom)

sample_kernel= outer( anno$findiv, anno$findiv, "==")

resdir=paste0(DATADIR,"sqtl_",normalization_approach,"_",permuted,"/")
dir.create(resdir)

checkpoint_dir=paste0(resdir,chrom,"_checkpoint/")
dir.create(checkpoint_dir)

results = map_interaction_qtl(input, genotype, geneloc, snploc, anno, sample_kernel, normalization_approach, permutation_approach, cisdist, checkpoint_dir) 

gz1 = gzfile(paste0(resdir,chrom,".txt.gz"),"w")
results %>% format(digits=5) %>% write.table(gz1, quote = F, row.names = F, col.names = T, sep="\t")
close(gz1)
