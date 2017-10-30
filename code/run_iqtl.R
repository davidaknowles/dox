
DATADIR="~/scailscratch/dox/"
library("dplyr")
library(data.table)
source("utils.R")
registerDoMC(16)

genotype=fread(paste0(DATADIR,"genotype.txt"))
setDF(genotype)
rownames(genotype)=genotype$V1
genotype$V1=NULL
colnames(genotype)=genotype[1,]
genotype=genotype[2:nrow(genotype),]

genotype=as.matrix(genotype)

sample_anno=read.table("../data/annotation.txt", header=T, stringsAsFactors = F)
errorCovariance=get_relatedness("../data/addSNP.coef.3671", unique(sample_anno$dbgap))

geneloc=read.table(paste0(DATADIR,"genelocGRCh38.txt"),header=T,stringsAsFactors = F)
snploc=read.table(paste0(DATADIR,"snploc.txt"),header=T,stringsAsFactors = F)

chrom=if (interactive()) "chr8" else commandArgs(trailingOnly = T)[1]

geneloc=geneloc[geneloc$chr==chrom,]
snploc=snploc[snploc$chr==chrom,]

stopifnot(all(as.character(snploc$snpid) %in% rownames(genotype) ))
genotype=genotype[as.character(snploc$snpid),]

input <- read.delim("../data/counts_log_cpm.txt.gz")

anno <- read.delim("../data/sample_annotation.txt", stringsAsFactors = F)

sample_anno=read.table("../data/annotation.txt", header=T, stringsAsFactors = F)

# mapping from cell-line ID to individual
dbgap=sample_anno$dbgap
names(dbgap)=sample_anno$cell_line
stopifnot(is.character(anno$individual))

colnames(input)=dbgap[anno$individual]

#input=remove_PCs(input, num_PCs_to_remove)
input=quantile_normalize(input)

anno$dbgap=as.character(dbgap[anno$individual])

require(rstan)
lmm=stan_model("lmm.stan")

genes=intersect(rownames(input),geneloc$geneid)
rownames(geneloc)=geneloc$geneid
cisdist=1e5
results=setNames( foreach(gene=genes, .errorhandling='pass') %do% {
  print(gene)
  y=input[gene,]
  cis_snps=snploc[ ((geneloc[gene,"left"]-cisdist) < snploc$pos) & ((geneloc[gene,"right"]+cisdist) > snploc$pos), "snpid" ]
  cis_snps=as.character(cis_snps)
  # cis_snp=as.character(cis_snps)[1]
  same_ind=outer(anno$dbgap, anno$dbgap, "==") * 1
  same_conc=outer(anno$conc, anno$conc, "==") * 1
  N=length(y)
  
  x_no_geno=list(diag(N),errorCovariance[ anno$dbgap, anno$dbgap ],same_ind,same_conc)
  data=list(N=N,x=x_no_geno,P=length(x_no_geno),y=y-mean(y))

  fit_no_geno=optimizing(lmm, data, as_vector=F)
  setNames( foreach(cis_snp=cis_snps, .errorhandling='pass') %dopar% {
    geno=genotype[cis_snp,anno$dbgap]
    #l=lm(y ~ geno + as.factor(anno$conc))
    #anno$geno=geno
    #lme(y ~ geno + as.factor(conc), anno, ~ 1|as.factor(dbgap), correlation = corSymm(, fixed=T))

    x_geno=c( x_no_geno, list(outer(geno,geno)) )
    data=list(N=N,x=x_geno,P=length(x_geno),y=y-mean(y))
    init=fit_no_geno$par
    init$s=c(init$s,0.01)
    fit_geno=optimizing(lmm, data, init=init, as_vector=F )
    
    interact=model.matrix(~geno:as.factor(conc),data=anno)
    interact=interact[,3:ncol(interact)]
    x_interact=c( x_geno, list(interact %*% t(interact) ) )
    data=list(N=N,x=x_interact,P=length(x_interact),y=y-mean(y))
    init=fit_geno$par
    init$s=c(init$s,0.01)
    fit_interact=optimizing(lmm, data, init=init, as_vector=F)
    
    V=foreach(p=seq_len(data$P-1), .combine = "+") %do% {fit_interact$par$s[p] * data$x[[p]]}
    xvx=t(interact) %*% solve(V, interact)
    beta_interact=solve( xvx, t(interact) %*% solve(V, y) )
    se_interact=sqrt(diag(solve( xvx, diag(nrow(xvx)) )))
    
    list(cis_snp=cis_snp, beta_interact=beta_interact,se_interact=se_interact,lrt=fit_interact$value - fit_geno$value)
  }, cis_snps )
}, genes )

save(results, file=paste0(DATADIR,"lmm_",chrom,".RData"))
