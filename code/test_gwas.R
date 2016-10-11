#num_PCs_to_remove=10
DATADIR="~/scailscratch/dox/"
library("dplyr")
library(data.table)
source("utils.R")
#registerDoMC(16)

genotype=fread(paste0(DATADIR,"genotype.txt"))
setDF(genotype)
rownames(genotype)=genotype$V1
genotype$V1=NULL
colnames(genotype)=genotype[1,]
genotype=genotype[2:nrow(genotype),]

genotype=as.matrix(genotype)

sample_anno=read.table("../data/annotation.txt", header=T, stringsAsFactors = F)
errorCovariance=get_relatedness("../data/addSNP.coef.3671", unique(sample_anno$findiv))

geneloc=read.table(paste0(DATADIR,"genelocGRCh38.txt"),header=T,stringsAsFactors = F)
snploc=read.table(paste0(DATADIR,"snploc.txt"),header=T,stringsAsFactors = F)

chrom=if (interactive()) "chr8" else commandArgs(trailingOnly = T)[1]

geneloc=geneloc[geneloc$chr==chrom,]
snploc=snploc[snploc$chr==chrom,]

input <- read.delim("../data/counts_log_cpm.txt.gz")

anno <- read.delim("../data/sample_annotation.txt", stringsAsFactors = F)

sample_anno=read.table("../data/annotation.txt", header=T, stringsAsFactors = F)

# mapping from cell-line ID to individual
findiv=sample_anno$findiv
names(findiv)=sample_anno$cell_line
stopifnot(is.character(anno$individual))

colnames(input)=findiv[anno$individual]

#input=remove_PCs(input, num_PCs_to_remove)
input=quantile_normalize(input)

anno$findiv=as.character(findiv[anno$individual])

require(rstan)
sm=stan_model("lmm.stan")

#gene="ENSG00000140396"
#cis_snp="6375587"

cis_snp="46760"
gene="ENSG00000049239"

y=input[gene,]

geno=genotype[cis_snp,anno$findiv]
print(ggplot(data.frame(y=y, geno=as.factor(geno), conc=anno$conc), aes(as.factor(conc), y, col=geno)) + geom_boxplot() + ggtitle(paste("Gene:",gene,"SNP:",cis_snp)) + ylab("Expression") + xlab("Dox concentration") + theme_bw(base_size=16))


geno=genotype[cis_snp,anno$findiv]

to_keep=!is.na(geno)

anno=anno[to_keep,]
geno=geno[to_keep]
y=y[to_keep]

same_ind=outer(anno$findiv, anno$findiv, "==") * 1
same_conc=outer(anno$conc, anno$conc, "==") * 1
N=length(y)

x_no_geno=list(diag(N),errorCovariance[ anno$findiv, anno$findiv ],same_ind,same_conc)
data=list(N=N,x=x_no_geno,P=length(x_no_geno),y=y-mean(y))

require(rstan)
sm=stan_model("lmm.stan")

fit_no_geno=optimizing(sm, data, as_vector=F)

x_geno=c( x_no_geno, list(outer(geno,geno)) )
data_geno=list(N=N,x=x_geno,P=length(x_geno),y=y-mean(y))
init=fit_no_geno$par
init$s=c(init$s,0.01)
fit_geno=optimizing(sm, data_geno, init=init, as_vector=F )

interact=model.matrix(~geno:as.factor(conc),data=anno)
interact=interact[,3:ncol(interact)]
x_interact=c( x_geno, list(interact %*% t(interact) ) )
data_interact=list(N=N,x=x_interact,P=length(x_interact),y=y-mean(y))
init=fit_geno$par
init$s=c(init$s,0.01)
fit_interact=optimizing(sm, data_interact, init=init, as_vector=F)

Lobs=fit_interact$value - fit_geno$value

V=foreach(i=seq_along(fit_geno$par$s), .combine = "+") %do% { fit_geno$par$s[i] * data_geno$x[[i]] }

L_bootstrap=foreach(i=1:1000, .combine = c) %dopar% {
  y_gen=MASS::mvrnorm(n = 1, numeric(N), V)
  
  x_geno=c( x_no_geno, list(outer(geno,geno)) )
  data_geno=list(N=N,x=x_geno,P=length(x_geno),y=y_gen)
  #init=fit_no_geno$par
  #init$s=c(init$s,0.01)
  init=fit_geno$par
  fit_geno_pb=optimizing(sm, data_geno, init=init, as_vector=F )
  
  interact=model.matrix(~geno:as.factor(conc),data=anno)
  interact=interact[,3:ncol(interact)]
  x_interact=c( x_geno, list(interact %*% t(interact) ) )
  data_interact=list(N=N,x=x_interact,P=length(x_interact),y=y_gen)
  init=fit_geno_pb$par
  init$s=c(init$s,0.01)
  fit_interact=optimizing(sm, data_interact, init=init, as_vector=F)
  fit_interact$value - fit_geno_pb$value
}

hist(L_bootstrap,100)
L_bootstrap_0=pmax(L_bootstrap,0)

qqplot( asinh(rchisq(1000,.5)), asinh( 2*L_bootstrap_0 ) )
abline(0,1)

pchisq(Lobs, 0.2, lower.tail = F)

lmm_with_fix=stan_model("lmm_with_fix.stan")

x_geno=c( x_no_geno, list(outer(geno,geno)) )
data_geno=list(N=N,x=x_geno,P=length(x_geno),y=y-mean(y))
init=fit_no_geno$par
init$s=c(init$s,0.01)
fit_geno=optimizing(sm, data_geno, init=init, as_vector=F )

interact=model.matrix(~geno:as.factor(conc),data=anno)
interact=interact[,3:ncol(interact)]
data_interact=list(N=N,x=x_geno,P=length(x_geno),xfix=interact,Pfix=ncol(interact),y=y-mean(y))
init=fit_geno$par
init$beta=numeric(data_interact$Pfix)
fit_interact=optimizing(lmm_with_fix, data_interact, init=init, as_vector=F)
lrt=2*(fit_interact$value - fit_geno$value)
pchisq(lrt, df=data_interact$Pfix, lower.tail = F)

L_fix_bootstrap=foreach(i=1:1000, .combine = c) %dopar% {
  y_gen=MASS::mvrnorm(n = 1, numeric(N), V)
  
  x_geno=c( x_no_geno, list(outer(geno,geno)) )
  data_geno=list(N=N,x=x_geno,P=length(x_geno),y=y_gen)
  #init=fit_no_geno$par
  #init$s=c(init$s,0.01)
  init=fit_geno$par
  fit_geno_pb=optimizing(sm, data_geno, init=init, as_vector=F )
  
  interact=model.matrix(~geno:as.factor(conc),data=anno)
  interact=interact[,3:ncol(interact)]
  data_interact=list(N=N,x=x_geno,P=length(x_geno),xfix=interact,Pfix=ncol(interact),y=y_gen)
  init=fit_geno_pb$par
  init$beta=numeric(data_interact$Pfix)
  fit_interact=optimizing(lmm_with_fix, data_interact, init=init, as_vector=F)
  
  fit_interact$value - fit_geno_pb$value
}

hist(L_fix_bootstrap,100)
qqplot( rchisq(1000,data_interact$Pfix), 2*L_fix_bootstrap)
abline(0,1)

