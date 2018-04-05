require(magrittr) 

require(glmnet)
library("dplyr")
library("tidyr")
library(data.table)
source("utils.R")

source("load_data.R")

geneloc=read.table(paste0(DATADIR,"genelocGRCh38.txt.gz"),header=T,stringsAsFactors = F)

cisdist=1e5
normalization_approach="log"

if (interactive()) {
  chrom="chr8"
} else {
  ca=commandArgs(trailingOnly = T)
  chrom=ca[1]
  registerDoMC(16)
}

geneloc=geneloc[geneloc$chr==chrom,]
snploc=snploc[snploc$chr==chrom,]

# if (interactive()) {
#    geneloc=head(geneloc,2)
# }

# stopifnot(all(as.character(snploc$snpid) %in% rownames(genotype) ))
# genotype=genotype[as.character(snploc$snpid),]

resdir=paste0(DATADIR,"predicted_expression_",normalization_approach,"/")
dir.create(resdir)

#' Simple SVD based imputation of missing genotypes
#' 
#' @param geno [samples] x [SNPs] genotype matrix (0/1/2)
#' @param prop_var Proportion of variance that the PCs should explain
#' 
#' @return Complete genotype matrix. 
easy_impute=function(geno, prop_var=0.95) {
  temp=geno
  temp=t(scale(t(geno)))
  temp[is.na(temp)]=0
  s=svd(temp)
  v=s$d^2/sum(s$d^2)
  to_use=cumsum(v)<prop_var
  s$d[!to_use]=0.0
  recon=s$u %*% fix_diag(s$d) %*% t(s$v)
  temp[is.na(geno)]=recon[is.na(geno)]
  temp=unscale(temp)
  stopifnot(max(abs(temp[!is.na(geno)]-geno[!is.na(geno)]))<1e-10)
  temp=round(temp)
  class(temp)="integer"
  temp
}

#' qqnorm without plotting
#' 
#' @param x Numeric vector input
#' @return Quantile normalized version
qqnorm_no_plot=function(x) {
  n=length(x)
  a=if (n<=10) (3.0/8.0) else 0.5
  qnorm( (rank(x)-a)/(n+1.0-2.0*a) )
}

#' Quantile normalize columns
quantile_normalize_cols=function(input) {
  apply(input, 2, qqnorm_no_plot)
}

#' Quantile normalize rows (to normal)
#' @import magrittr
quantile_normalize=function(input) {
  apply(input, 1, qqnorm_no_plot) %>% t()
}

if (normalization_approach=="qq") {
  input=quantile_normalize(input)
}

genes=intersect(rownames(input),geneloc$geneid)
rownames(geneloc)=geneloc$geneid

anno = anno %>% mutate(condition=as.factor(conc))

res=foreach(gene=genes) %dopar% {
  cis_snps=snploc %>% 
    filter(chr==geneloc[gene,"chr"], 
           (geneloc[gene,"left"]-cisdist) < pos,  
           (geneloc[gene,"right"]+cisdist) > pos) %>%
    .$snpid %>%
    as.character()
  cat(gene,length(cis_snps)," cis snps\n")
  
  cis_snps=intersect(rownames(genotype), cis_snps)
  
  if (length(cis_snps)==0) return(NULL)
  
  y=input[gene,] %>% as.numeric
  # y=y-mean(y)
  
  imp_geno=easy_impute(genotype[cis_snps,,drop=F])
  
  x=imp_geno[,anno$dbgap,drop=F] %>% t() %>% as.data.frame() %>% cbind(conc=factor(anno$conc))
  x=model.matrix( ~ .*conc, data=x )
  
  cv=cv.glmnet(x, y, alpha=0.5, keep=T)
  cv$fit.preval[,which.min(cv$cvm)]
} %>% set_names(genes) 

res=res[!sapply(res,is.null)]

gz=gzfile(paste0(resdir,chrom,".txt.gz"),"w")
res %>% as.data.frame() %>% write.table(gz, quote=F, col.names=T, sep="\t")
close(gz)

