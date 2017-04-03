require(leafcutter)
library(abind)
require(rstan)
require(data.table)
require(doMC)

source("utils.R")

require(reshape2)
registerDoMC(detectCores()-1)

counts=read.table("../data/counts.txt.gz", check.names = F)

anno <- read.delim("../data/sample_annotation.txt", stringsAsFactors = F)

library_size_model=stan_model("~/Dropbox/eagle/eagle/beta_binomial_models/library_size.stan")
counts=t(counts)
dat=list(N=nrow(counts), P=ncol(counts), gene_counts=counts, concShape=1.0001, concRate=0.0001)

gene_means=colMeans(counts)
gene_means=gene_means/sum(gene_means)
library_size_guess=rowMeans(sweep(counts, 2, gene_means, "/"))

mean( (counts - outer(library_size_guess, gene_means))^2 ) / mean(counts^2)

init=list( library_size=library_size_guess, gene_means=gene_means)
o=optimizing(library_size_model, dat, init=init, as_vector=F)

library_size_fit=o$par$library_size

plot(library_size_guess, library_size_fit); abline(0,1)

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

missing_from_as=setdiff(rownames(counts), colnames(n))
length(missing_from_as) # 25. No idea what happened to these. 
mean(counts[missing_from_as,], na.rm=T) 
mean(counts[intersect(rownames(counts), colnames(n)),], na.rm=T)

missing_from_total=setdiff(colnames(n),rownames(counts))
length(missing_from_total) # 12
mean(n[,intersect(rownames(counts), colnames(n))]) # ~=10. 
mean(n[,missing_from_total]) # ~=1. These have very small library size. 

tokeep=intersect(rownames(counts), colnames(n))

rownames(samp_meta)=colnames(n)
alt=alt[,tokeep]
n=n[,tokeep]
samp_meta=samp_meta[tokeep,]

is_het=stan_model("~/Dropbox/eagle/eagle/beta_binomial_models/is_het.stan")

cast_counts=function(g) {
  temp=cbind(anno, val=g)
  d=dcast(temp, individual ~ conc  , value.var="val")
  rownames(d)=d$individual
  d$individual=NULL
  d
}

library_size_mat=cast_counts(library_size_fit / mean(library_size_fit))

class(missing)="integer"

library_size_mat[is.na(library_size_mat)]=0

filter_data=function(gene) {
  
  snps=gene_snps[gene_snps$gene==gene,"snp"]
  i=snpmeta %in% snps
  
  a=alt[i,,drop=F]
  a=do.call(abind,c( apply(a, 1, function(g) {
    temp=cbind(samp_meta, a=g)
    d=dcast(temp, samp ~ timepoint , value.var="a")
    rownames(d)=d$samp
    d$samp=NULL
    d
  }), along=3 ))
  a[is.na(a)]=0
  
  nh=n[i,,drop=F]
  nh=do.call(abind,c( apply(nh, 1, function(g) {
    temp=cbind(samp_meta, n=g)
    d=dcast(temp, samp ~ timepoint , value.var="n")
    rownames(d)=d$samp
    d$samp=NULL
    d
  }), along=3 ))
  nh[is.na(nh)]=0 # need this to account for missing samples
  
  dummy=foreach(snp_index=seq_len(dim(nh)[3])) %do% {
  #foreach(snp_index=3:8) %do% {
    as=a[,,snp_index]
    nhh=nh[,,snp_index]
    ind_to_keep=rowSums(nhh)>0
    treat_to_keep=colSums(nhh)>0
    
    as=as[ind_to_keep,treat_to_keep,drop=F]
    nhh=nhh[ind_to_keep,treat_to_keep,drop=F]
    
    o=  optimizing(is_het, dat=list(N=nrow(nhh), T=ncol(nhh), errorRate=0.01, concShape=1.001, concRate=0.001, ys=as, ns=nhh), as_vector=F)
    eo=exp(o$par$probs)
    pr=sweep(eo, 1, rowSums(eo), "/")
    
    homo=pr[,1]<0.95
    
    cat("Removing",sum(homo),"individual(s) from SNP",dimnames(a)[[3]][snp_index],"\n")
    
    nh[which(ind_to_keep)[homo],treat_to_keep,snp_index]=0
    a[which(ind_to_keep)[homo],treat_to_keep,snp_index]=0
  }
  
  snp_to_keep=apply(nh>0, 3, any)
  a=a[,,snp_to_keep,drop=F]
  nh=nh[,,snp_to_keep,drop=F]
  
  counts_here=cast_counts(counts[,gene])
  counts_here[is.na(counts_here)]=0
  
  list(ys=a,ns=nh,gene_counts=counts_here)
}


genotype=fread(paste0("zcat < ../data/genotype.txt.gz"))
setDF(genotype)
rownames(genotype)=genotype$V1
genotype$V1=NULL
colnames(genotype)=genotype[1,]
genotype=genotype[2:nrow(genotype),]
genotype=as.matrix(genotype)

sample_anno=read.table("../data/annotation.txt", header=T, stringsAsFactors = F)
findiv=sample_anno$findiv
names(findiv)=sample_anno$cell_line
class(findiv)="character"
stopifnot(is.character(anno$individual))

genes=intersect( unique(gene_snps$gene), colnames(counts) )

geneloc=read.table("../data/genelocGRCh38.txt.gz",header=T,stringsAsFactors = F)
snploc=read.table("../data/snploc.txt.gz",header=T,stringsAsFactors = F)

snploc$chrpos=with(snploc, paste(chr,pos,sep=":"))

genes=intersect(genes,geneloc$geneid)
rownames(geneloc)=geneloc$geneid

eagle2=stan_model("~/Dropbox/eagle/eagle/beta_binomial_models/eagle2.stan")

gene="ENSG00000166435"
cis_snp="8198371"

gene="ENSG00000163682"
cis_snp="3039070"

gene="ENSG00000066827"
cis_snp="6628114"

gene="ENSG00000242110"
cis_snp="3883582"

gene="ENSG00000136371"
cis_snp="10170248"

gene="ENSG00000123178"
cis_snp="9216771"

lrt_res=read.table("../data/lrt-summary.txt.gz", stringsAsFactors = F, header=T)
top_iqtl= order(lrt_res$p)[3] # 3, 
gene=lrt_res$gene[top_iqtl]
cis_snp=as.character(lrt_res$snp[top_iqtl])
lrt_res[top_iqtl,]

cisdist=1e5
errorhandling=if (interactive()) 'stop' else 'pass'
results=setNames( foreach(gene=genes, .errorhandling=errorhandling) %do% {
  print(gene)
  
  dat=filter_data(gene)
  
  stopifnot(all(rownames(dat$gene_counts)==dimnames(dat$ns)[[1]]))
  
  cis_snps=snploc[ (geneloc[gene,"chr"]==snploc$chr) & ((geneloc[gene,"left"]-cisdist) < snploc$pos) & ((geneloc[gene,"right"]+cisdist) > snploc$pos), "snpid" ]
  cis_snps=as.character(cis_snps)
  
  imp_geno=easy_impute(genotype[cis_snps,])
  
  res=rbindlist( foreach(cis_snp=cis_snps, .errorhandling=errorhandling) %dopar% {
    geno=imp_geno[cis_snp,findiv[rownames(dat$gene_counts)] ]
    
    # matrix[N,P] x_1[T]. 
    dat$C=ncol(dat$gene_counts) # num conditions
    dat$N=nrow(dat$gene_counts)
    dat$P=dat$C+1 # mean for each group + haplo(geno)type
    temp=matrix(0,dat$N,dat$C)
    temp[,1]=0.5
    dat$x_1=foreach(condition=1:dat$C) %do% {
      tt=temp
      tt[,condition]=0.5
      cbind( tt, geno>=1 ) # arbitrary hap assignment
    }
    dat$x_2=foreach(condition=1:dat$C) %do% {
      tt=temp
      tt[,condition]=0.5
      cbind( tt, geno>=2 )
    }
    
    data_geno=c(dat, list( library_size=library_size_mat, N=length(geno), K=dim(dat$ys)[3], T=dat$C, P=dat$P, concShape=1.001, concRate=0.001) )
    #data_geno$library_size[ (dat$gene_counts / library_size_mat) > 2000 ]=0
    
    melt_it=function(g, value.name){
      g$ind=rownames(g)
      melt(g, id.vars = "ind", variable.name = "treat", value.name=value.name)
    }
    
    lin_data=melt_it(dat$gene_counts,"raw")
    lin_data$lib_size=melt_it(library_size_mat, "value")$value
    lin_data$count=with(lin_data, raw/lib_size)
    lin_data$geno=geno[ findiv[ lin_data$ind ] ]
    lin_data=lin_data[!is.na(lin_data$count),]
    lin_fit=lm(count ~ treat + geno, data=lin_data)
    lin_data$fitted=fitted(lin_fit)
    
    get_conc_init=function(fit) {
      fitt=fitted(fit)
      fitt[fitt <= 1]=1
      1/mean(1/with(lin_data, fitt * lib_size / ((1+raw) * (residuals(fit)^2 * lib_size / fitt - 1))))
    }
   # ggplot(lin_data, aes(treat, col=as.factor(geno))) + geom_boxplot(aes(y=fitted), outlier.shape = NA) + geom_point(aes(y=count), position = position_dodge(width=2/3),alpha=.5)
    
    #init=list(beta=c(1000,numeric(dat$P-1)), conc=1, prior=array(0.5+numeric(data_geno$K)))
    init=list(beta=coef(lin_fit), conc=get_conc_init(lin_fit), prior=array(0.5+numeric(data_geno$K)))
    if (any(is.na(init$beta))) return(NULL)
    
    fit_genos=foreach(sc=c(0.1,1,10)) %do% {
      init_temp=init
      init_temp$conc=init$conc * sc
      optimizing(eagle2, data_geno, init=init_temp, verbose=T, as_vector=F, algorithm=algorithm )
    }
    fit_geno=fit_genos[[ which.max(sapply(fit_genos, function(g) g$value)) ]]
    # model.matrix( ~cond*hap, data=data.frame( cond=as.factor(rep(1:5,10)), hap=rep(0:1,25)) )
    
    data_interact=data_geno
    data_interact$P=dat$C+1+(dat$C-1)
    data_interact$x_1=foreach(condition=1:dat$C) %do% {
      tt=temp
      tt[,condition]=.5
      hap=geno>=1
      cbind( tt, hap, sweep(tt[,2:dat$C],1,hap,"*") ) # arbitrary hap assignment
    }
    data_interact$x_2=foreach(condition=1:dat$C) %do% {
      tt=temp
      tt[,condition]=.5
      hap=geno>=2
      cbind( tt, hap, sweep(tt[,2:dat$C],1,hap,"*") )
    }
    df=data_interact$P-data_geno$P
    
    lin_fit_interact=lm(count ~ treat * geno, data=lin_data)
    
    #lin_data$interact_fitted=fitted(lin_fit_interact)
    #ggplot(lin_data, aes(treat, interact_fitted, col=as.factor(geno))) + geom_boxplot() + ylim(0,3500)
    
    #init=fit_geno$par
    #init$beta=c(init$beta,numeric(df))
    #init$beta=coef(lin_fit_interact)
    init_interact=list(beta=coef(lin_fit_interact), conc=get_conc_init(lin_fit_interact), prior=array(0.5+numeric(data_geno$K)))
    if (any(is.na(init_interact$beta))) return(NULL)
    
    fit_interacts=foreach(sc=c(0.1,1,10)) %do% {
      init_temp=init_interact
      init_temp$conc=init_interact$conc * sc
      optimizing(eagle2, data_interact, init=init_temp, verbose=T, as_vector=F, algorithm=algorithm )
    }
    fit_interact=fit_interacts[[ which.max(sapply(fit_interacts, function(g) g$value)) ]]
    
    lrt=2.0*(fit_interact$value - fit_geno$value)
    
    lin_data$ynorm=qqnorm(lin_data$count, plot.it = F)$x
    p_qnorm=anova(lm(ynorm ~ treat + geno, data=lin_data),lm(ynorm ~ treat * geno, data=lin_data))[["Pr(>F)"]][2]
    #ggplot(lin_data, aes(treat, ynorm, col=as.factor(geno))) + geom_boxplot() + geom_point(position = position_dodge(width=.6))
    #ggplot(lin_data, aes(treat, count, col=as.factor(geno))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_dodge(width=4/5),alpha=.5) + theme_bw()
    lin_data$interact_fitted=fitted(lin_fit_interact)
    #ggplot(lin_data, aes(treat, interact_fitted, col=as.factor(geno))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_dodge(width=.6))
    data.frame(gene=gene, cis_snp=cis_snp, lrt=lrt, df=df, nln10p=-pchisq(lrt, df, lower.tail = F, log.p = T)/log(10), plmm=lrt_res$p[top_iqtl],  p_lm=anova(lin_fit_interact, lin_fit)[["Pr(>F)"]][2], p_qnorm=p_qnorm, init_conc=get_conc_init(lin_fit), fit_conc=fit_geno$par$conc)
  } )
}, genes )


df=as.data.frame(dat$gene_counts/library_size_mat)
df$ind=rownames(df)
m=melt(df,id.vars = "ind",variable.name = "treat", value.name = "count")
m$geno=as.factor(geno[ findiv[ m$ind ] ])
m$treat=as.factor(m$treat)
#ggplot(m, aes(interaction(geno,treat), count)) + geom_boxplot()
ggplot(m, aes(treat, count, col=geno)) + geom_boxplot()

ar=dat$ys / dat$ns

snp_coverage=apply(dat$ns, 3, sum)
which_snp=which.max(snp_coverage)
#which_snp=2
ar=as.data.frame(ar[ ,, which_snp ])

ar$ind=rownames(ar)

m=melt(ar, id.vars = "ind", variable.name = "treat", na.rm = T)
m$geno=as.factor(geno[ findiv[ m$ind ] ])



ggplot(m, aes(treat, value, col=geno)) +geom_point(position = position_jitter(0.2))

gzf=gzfile("../data/eagle_joint.txt.gz","w")
write.table(data.frame(gene=genes,p=p), file=gzf, row.names = F, quote = F, sep="\t")
close(gzf)
cat("Done!\n")
stop()

require(ggplot2)

nhhh=dat$ns
temp=kronecker(library_size_mat==0, array(T, dim=c(1,1,dim(nhhh)[3])))
class(temp)="logical"
nhhh[temp]=NA

nhh=as.data.frame(nhhh[,,1])
df=as.data.frame(nhh)
df$ind=rownames(df)
m=melt(df, id.vars = "ind", variable.name = "treat", na.rm = T)
m$geno=geno[ findiv[m$ind] ]
m=m[m$value>0,]

ggplot(m, aes(treat, value, col=as.factor(geno))) + geom_boxplot()
