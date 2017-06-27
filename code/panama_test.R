
DATADIR="~/scailscratch/dox/"
library("dplyr")
library("tidyr")
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

findiv[ findiv==160001 ]=106411

anno$findiv=as.character(findiv[anno$individual])

require(rstan)
panama_test=stan_model("panama_test.stan")

genes=intersect(rownames(input),geneloc$geneid)
rownames(geneloc)=geneloc$geneid
cisdist=1e5
errorhandling=if (interactive()) 'stop' else 'remove'

same_ind=outer(anno$findiv, anno$findiv, "==") * 1
same_conc=outer(anno$conc, anno$conc, "==") * 1

eig_K=eigen(readRDS("../data/Kern.rds"))

results=foreach(gene=genes, .errorhandling=errorhandling, .combine = bind_rows) %do% {
  
  print(gene)
  y=input[gene,]
  cis_snps=snploc[ ((geneloc[gene,"left"]-cisdist) < snploc$pos) & ((geneloc[gene,"right"]+cisdist) > snploc$pos), "snpid" ]
  cis_snps=as.character(cis_snps)

  imp_geno=easy_impute(genotype[cis_snps,])
  # cis_snp=as.character(cis_snps)[1]
  
  N=length(y)
  
  intercept_only=matrix(1,ncol=1,nrow=N)
  data=list(N=N,x=intercept_only,P=1,y=y-mean(y), U_transpose=t(eig_K$vectors), lambda=eig_K$values)

  init=list(sigma2=0.1, sigma2_k=1.0, beta=array(0.0))

  fit_no_geno=optimizing(panama_test, data, init=init, as_vector=F)
  foreach(cis_snp=cis_snps, .errorhandling=errorhandling, .combine = bind_rows) %dopar% {
    geno=imp_geno[cis_snp,anno$findiv]
    if (sum(imp_geno[cis_snp,]) < 5.0) return(NULL)
    #l=lm(y ~ geno + as.factor(anno$conc))
    #anno$geno=geno
    #lme(y ~ geno + as.factor(conc), anno, ~ 1|as.factor(findiv), correlation = corSymm(, fixed=T))

    data$x=cbind( intercept_only, geno )
    data$P=ncol(data$x)
    init=fit_no_geno$par
    init$beta=c(init$beta,0.0)
    
    fit_geno=optimizing(panama_test, data, init=init, as_vector=F )
    
    interact=model.matrix(~geno:as.factor(conc),data=anno)
    interact=interact[,3:ncol(interact)]
    data_interact=data
    data_interact$x=cbind( intercept_only, geno, interact )
    data_interact$P=ncol(data_interact$x)
  
    init=fit_geno$par
    init$beta=c(init$beta,numeric(ncol(interact)))
    fit_interact=optimizing(panama_test, data_interact, init=init, as_vector=F)
    
    #lrt_interact=2.0*(fit_interact$value - fit_geno$value)
    #df_interact=ncol(interact)
    # pchisq(lrt,df,lower.tail=F)
    
    data.frame(gene=gene, cis_snp=cis_snp, l0=fit_no_geno$value, l_geno=fit_geno$value, l_interact=fit_interact$value, df=ncol(interact)  )
  }
  
}

resdir=paste0("~/dagscratch/dox/panama_",normalization_approach,"/")
dir.create(resdir)

gz1 = gzfile(paste0(resdir,chrom,".txt.gz"),"w")
results %>% format(digits=5) %>% write.table(gz1, quote = F, row.names = F, col.names = T, sep="\t")
close(gz1)
