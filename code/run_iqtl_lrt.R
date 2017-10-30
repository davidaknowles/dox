
DATADIR="~/scailscratch/dox/"
library("dplyr")
library(data.table)
source("utils.R")
registerDoMC(16)

genotype=fread("zcat < ../data/genotype.txt.gz", data.table = F, header = T)

rownames(genotype)=genotype$snpid
genotype$snpid=NULL

genotype=as.matrix(genotype)

sample_anno=read.table("../data/annotation.txt", header=T, stringsAsFactors = F)
if (F) {
    errorCovariance=get_relatedness("../data/addSNP.coef.3671", unique(sample_anno$findiv))
    saveRDS( errorCovariance, file="../data/error_covariance.Rds" )
} else { errorCovariance = readRDS( "../data/error_covariance.Rds" ) }

geneloc=read.table(paste0(DATADIR,"genelocGRCh38.txt"),header=T,stringsAsFactors = F)
snploc=read.table(paste0(DATADIR,"snploc.txt"),header=T,stringsAsFactors = F)

if (interactive()) {
  chrom="chr8"
  normalization_approach="qq"
} else {
  ca=commandArgs(trailingOnly = T)
  chrom=ca[2]
  normalization_approach=ca[1]
}

geneloc=geneloc[geneloc$chr==chrom,]
snploc=snploc[snploc$chr==chrom,]

stopifnot(all(as.character(snploc$snpid) %in% rownames(genotype) ))
genotype=genotype[as.character(snploc$snpid),]

input <- read.delim("../data/counts_log_cpm.txt.gz")

anno <- read.delim("../data/sample_annotation.txt", stringsAsFactors = F)

sample_anno=read.table("../data/annotation.txt", header=T, stringsAsFactors = F)

# mapping from cell-line ID to individual
findiv=sample_anno$findiv
names(findiv)=sample_anno$cell_line
stopifnot(is.character(anno$individual))

colnames(input)=findiv[anno$individual]

#input=remove_PCs(input, num_PCs_to_remove)
if (normalization_approach=="qq") {
  input=quantile_normalize(input)
} else if (normalization_approach=="log") {
  input=input %>% t %>% scale %>% t
} else if (normalization_approach=="linear") {
  input=(2^input) %>% t %>% scale %>% t
}

findiv[ findiv=="7440_4ce2" ]="3e07_41cd"

anno$findiv=as.character(findiv[anno$individual])

require(rstan)
lmm=stan_model("lmm.stan")
lmm_with_fix=stan_model("lmm_with_fix.stan")

genes=intersect(rownames(input),geneloc$geneid)
rownames(geneloc)=geneloc$geneid
cisdist=1e5
errorhandling=if (interactive()) 'stop' else 'pass'

results=setNames( foreach(gene=genes, .errorhandling=errorhandling) %do% {
  print(gene)
  y=input[gene,]
  cis_snps=snploc[ ((geneloc[gene,"left"]-cisdist) < snploc$pos) & ((geneloc[gene,"right"]+cisdist) > snploc$pos), "snpid" ]
  cis_snps=as.character(cis_snps)

  imp_geno=easy_impute(genotype[cis_snps,])
  # cis_snp=as.character(cis_snps)[1]
  same_ind=outer(anno$findiv, anno$findiv, "==") * 1
  same_conc=outer(anno$conc, anno$conc, "==") * 1
  N=length(y)
  
  x_no_geno=list(diag(N),errorCovariance[ anno$findiv, anno$findiv ],same_ind,same_conc)
  data=list(N=N,x=x_no_geno,P=length(x_no_geno),y=y-mean(y))

  fit_no_geno=optimizing(lmm, data, as_vector=F)
  setNames( foreach(cis_snp=cis_snps, .errorhandling=errorhandling) %dopar% {
    geno=imp_geno[cis_snp,anno$findiv]
    #l=lm(y ~ geno + as.factor(anno$conc))
    #anno$geno=geno
    #lme(y ~ geno + as.factor(conc), anno, ~ 1|as.factor(findiv), correlation = corSymm(, fixed=T))

    x_geno=c( x_no_geno, list(outer(geno,geno)) )
    data=list(N=N,x=x_geno,P=length(x_geno),y=y-mean(y))
    init=fit_no_geno$par
    init$s=c(init$s,0.01)
    fit_geno=optimizing(lmm, data, init=init, as_vector=F )
    
    interact=model.matrix(~geno:as.factor(conc),data=anno)
    interact=interact[,3:ncol(interact)]
    data_interact=list(N=N,x=x_geno,P=length(x_geno),xfix=interact,Pfix=ncol(interact),y=y-mean(y))
    init=fit_geno$par
    init$beta=numeric(data_interact$Pfix)
    fit_interact=optimizing(lmm_with_fix, data_interact, init=init, as_vector=F)
    
    lrt=2.0*(fit_interact$value - fit_geno$value)
    df=data_interact$Pfix

    list(cis_snp=cis_snp, beta_interact=fit_interact$par$beta, lrt=lrt, df=df, p=pchisq(lrt, df, lower.tail = F))
  }, cis_snps )
}, genes )

resdir=paste0("~/dagscratch/dox/lrt_imp_",normalization_approach,"/")
dir.create(resdir)
save(results, file=paste0(resdir,chrom,".RData"))
