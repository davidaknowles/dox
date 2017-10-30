require(leafcutter)
library(abind)
require(rstan)
require(data.table)
require(doMC)

source("utils.R")

require(reshape2)
registerDoMC(detectCores()-1)

anno <- read.delim("../data/sample_annotation.txt", stringsAsFactors = F)

gene_snps=fread("zcat < ../data/gene_snp_mapping.txt.gz")
setDF(gene_snps)

dat=fread("zcat < ../data/ase.txt.gz")
setDF(dat)

dat$snp=with(dat, paste(chr,pos,sep=":"))

ref=reshape2::dcast(dat, snp ~ sample, value.var = "r")
rownames(ref)=ref$snp
ref$snp=NULL
ref=as.matrix(ref)

alt=reshape2::dcast(dat, snp ~ sample, value.var = "y")
rownames(alt)=alt$snp
alt$snp=NULL
alt=as.matrix(alt)

stopifnot(all(colnames(alt)==colnames(ref)))
stopifnot(all(rownames(alt)==rownames(ref)))

alt[is.na(alt)]=0
ref[is.na(ref)]=0

n=alt+ref

samp_meta=setNames( as.data.frame(do.call(rbind,strsplit(colnames(alt),"_"))), c("samp","timepoint") )

snpmeta=rownames(alt)

colnames(n)=with(samp_meta, paste(samp,timepoint,sep=":"))
colnames(alt)=colnames(n)

rownames(samp_meta)=colnames(n)


genotype=fread(paste0("zcat < ../data/genotype.txt.gz"))
setDF(genotype)
rownames(genotype)=genotype$V1
genotype$V1=NULL
colnames(genotype)=genotype[1,]
genotype=genotype[2:nrow(genotype),]
genotype=as.matrix(genotype)

sample_anno=read.table("../data/annotation.txt", header=T, stringsAsFactors = F)
dbgap=sample_anno$dbgap
names(dbgap)=sample_anno$cell_line
class(dbgap)="character"
stopifnot(is.character(anno$individual))

genes=intersect( unique(gene_snps$gene), colnames(counts) )

geneloc=read.table("../data/genelocGRCh38.txt.gz",header=T,stringsAsFactors = F)
snploc=read.table("../data/snploc.txt.gz",header=T,stringsAsFactors = F)

snploc$chrpos=with(snploc, paste(chr,pos,sep=":"))


is_het=stan_model("~/Dropbox/eagle/eagle/beta_binomial_models/is_het.stan")

cast_counts=function(g) {
  temp=cbind(samp_meta, a=g)
  d=dcast(temp, samp ~ timepoint , value.var="a")
  rownames(d)=d$samp
  d$samp=NULL
  d
}

library_size_mat[is.na(library_size_mat)]=0

genotyped=intersect(snpmeta, snploc$chrpos)
alt=alt[genotyped,]
n=n[genotyped,]

res=rbindlist( foreach(snp_index=seq_along(genotyped[1:500])) %dopar% {
  
  as_full=cast_counts(alt[snp_index,])
  nh_full=cast_counts(n[snp_index,])
  
  to_keep=rowSums(nh_full,na.rm = T) > 3
  
  if (to_keep==0) return(NULL)
  
  as=as_full[to_keep,,drop=F]
  nh=nh_full[to_keep,,drop=F]
  
  as[is.na(as)]=0
  nh[is.na(nh)]=0
  
  o=optimizing(is_het, dat=list(N=nrow(nh), T=ncol(nh), errorRate=0.01, concShape=1.001, concRate=0.001, ys=as, ns=nh), as_vector=F)
  eo=exp(o$par$probs)
  pr=sweep(eo, 1, rowSums(eo), "/")
  
  het=logical(nrow(as_full))
  het[to_keep]=pr[,1]>0.95
  data.frame(  het=het, geno=genotype[ rownames(alt)[snp_index] , dbgap[ rownames(as_full) ] ] )
} )

setDF(res)

res$geno_het=res$geno==1

ta=table(res[,c("het","geno_het")])
fpr=ta["TRUE","FALSE"] / sum(ta["TRUE",])
fnr=ta["FALSE","TRUE"] / sum(ta[,"TRUE"])
