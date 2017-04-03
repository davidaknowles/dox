require(leafcutter)
library(abind)
require(data.table)
require(doMC)
registerDoMC(detectCores()-1)

source("~/Dropbox/eagle/eagle/beta_binomial_models/bb_glm_flips_rep_prior_gene.R")

dat=fread("zcat < ../data/ase.txt.gz")
setDF(dat)

gene_snps=fread("zcat < ../data/gene_snp_mapping.txt.gz")
setDF(gene_snps)

dat$snp=with(dat, paste(chr,pos,sep=":"))

require(reshape2)

ref=dcast(dat, snp ~ sample, value.var = "r")
rownames(ref)=ref$snp
ref$snp=NULL
ref=as.matrix(ref)

alt=dcast(dat, snp ~ sample, value.var = "y")
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
  
  rat=a/nh
  rat[nh==0]=0
  invalid=(nh<4) | (rat< 0.01) | (rat> 0.99)
  invalid_geno=apply(invalid, c(1,3), any)
  ind_to_keep=apply(!invalid_geno, 1, any)
  snp_to_keep=apply(!invalid_geno, 2, any)
  if (sum(ind_to_keep)<2 | sum(snp_to_keep)==0) return(NULL)
  a=a[ind_to_keep,,snp_to_keep,drop=F]
  nh=nh[ind_to_keep,,snp_to_keep,drop=F]
  invalid_geno=invalid_geno[ind_to_keep,snp_to_keep,drop=F]
  to_0=aperm( abind( foreach(ggg=seq_len(dim(a)[2])) %do% invalid_geno, along = 3 ), c(1,3,2) )
  
  a[to_0]=0
  nh[to_0]=0
  nh[is.na(nh)]=0
  
  list(a=a,nh=nh)
}

genes=unique(gene_snps$gene)
p=foreach(gene=genes, .combine=c) %dopar% {
  
  filtered_data=filter_data(gene)
  if (is.null(filtered_data)) return(NA)
  
  #tryCatch({ 
     res <- betabinomial_glm_flips_rep_gene(filtered_data$a,filtered_data$nh,verbose=F, iterations = 5000) 
      res$lrtp[1]
 #     } , error=function(g) NA )
}

gzf=gzfile("../data/eagle_gene_1conc.txt.gz","w")
write.table(data.frame(gene=genes,p=p), file=gzf, row.names = F, quote = F, sep="\t")
close(gzf)
cat("Done!\n")
stop()

require(ggplot2)

#res=read.table("../data/eagle_gene.txt.gz",sep="\t",header=T, stringsAsFactors = F)

res=read.table("../data/eagle_gene_1conc.txt.gz",sep="\t",header=T, stringsAsFactors = F)


rownames(snpmeta)=snpmeta$variantID

require(dplyr)
gene_id=gene_snps %>% group_by(gene) %>% summarize(id=do.call(paste,as.list(snp)))
class(gene_id)="data.frame"
back_track=gene_id %>% group_by(id) %>% summarise(gene=do.call(paste,as.list(gene)))
class(back_track)="data.frame"
rownames(back_track)=back_track$id
gene_id$ambiguous=back_track[gene_id$id, "gene"]

unique_gene=gene_id[!duplicated(gene_id$id),]
rownames(res)=res$gene
res=res[unique_gene$gene,]
res$ambiguous=unique_gene$ambiguous

res$q=bh(res$p)
nsig=sum(res$q<.1, na.rm=T)
multiqq(list(gene=res$p))



pdf("../figures/eagle_gene1.pdf",width=12,height=10)
foreach(i=order(res$p)[1:nsig]) %do% {
  gene=res$gene[i]
  dat=filter_data(gene)
  melted=melt(dat$a/dat$nh)
  melted_n=melt(dat$nh)
  tokeep=!is.na(melted$value)
  melted=cbind(melted[tokeep,], n=melted_n[tokeep,"value"])
  colnames(melted)=c("Individual","TimePoint","SNP","AllelicRatio","Coverage")
  melted$TimePoint=as.numeric(melted$TimePoint)
  
  the_title=paste0(res$ambiguous[i]," p=",format(res$p[i],digits=1)," q=",format(res$q[i],digits=1))
  
  #levels(melted$SNP)=snpmeta[levels(melted$SNP),"position"]
  melted$inter=interaction( melted$TimePoint, melted$SNP )
  melted$ii=as.numeric(melted$inter)
  print( ggplot(melted, aes( inter, AllelicRatio, label=Individual, col=SNP, shape=Individual)) + geom_point( aes( size=Coverage)) + ylim(0,1) + theme_bw(base_size = 16) + geom_line(aes(ii,AllelicRatio),alpha=.5)  + scale_shape_manual(values=seq_along(levels(melted$Individual))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Treatment.SNP") + ggtitle(the_title) ) 
  
  NULL
}
dev.off()