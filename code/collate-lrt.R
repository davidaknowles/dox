require(ggplot2)
require(doMC)

DATADIR="~/scailscratch/dox/"
library("dplyr")
library(data.table)
source("utils.R")

genotype=fread("zcat < ../data/genotype.txt.gz", data.table=F)

rownames(genotype)=genotype$V1
genotype$V1=NULL
colnames(genotype)=genotype[1,]
genotype=genotype[2:nrow(genotype),]

genotype=as.matrix(genotype)

sample_anno=read.table("../data/annotation.txt", header=T, stringsAsFactors = F)

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

chroms=c(paste0("chr",1:22),"chrX")

pdf("../figures/panama_hits.pdf",width=7,height=5)
#pdf("../figures/lrt-iqtl-imp-new.pdf",width=7,height=5)
#threshold=0.0005
#threshold=0.00189
threshold=0.0
res=foreach(chrom=chroms, .combine = bind_rows) %do% {
  print(chrom)
  #fn=paste0("~/dagscratch/dox/lrt_imp/lrt_imp_",chrom,".RData")
  #fn=paste0("~/scailscratch/dox/lrt_imp_log/",chrom,".RData")
  fn=paste0("~/scailscratch/dox/panama_qq/",chrom,".RData")
  if (!file.exists(fn)) {
      cat("No file:",fn,"\n")
      return(NULL)
  }
  load(fn)
  foreach(gene=names(results), .combine=bind_rows) %dopar% {
    result=results[[gene]]
    if ("error" %in% class(result)) return(NULL)
    ps=foreach(cis_snp=names(result), .combine=c) %do% {
        res=result[[cis_snp]]
        if (!("error" %in% class(res))) res$p else NA
    }
    cis_snp=names(result)[which.min(ps)]
    bfp=min(ps,na.rm=T) * sum(!is.na(ps))
    if (!is.na(bfp)) if (bfp < threshold) {
      print(gene)
      y=input[gene,]
      geno=genotype[cis_snp,anno$findiv]
      print(data.frame(y=y, geno=as.factor(geno), conc=anno$conc) %>% filter(!is.na(geno)) %>% ggplot(aes(as.factor(conc), y, col=geno)) + geom_boxplot() + ggtitle(paste("Gene:",gene,"SNP:",cis_snp)) + ylab("Expression") + xlab("Dox concentration") + theme_bw(base_size=16))
      print(c(cis_snp, gene))
    }
    data.frame( gene=gene, snp=cis_snp, p=bfp )
  } 
  
} 
dev.off() 

#q=p.adjust(res, method="BH")
#sum(q<.05,na.rm = T)
# q<0.05 ~ p<0.0005

gzf=gzfile("../data/panama-summary.txt.gz","w")
write.table(res, file=gzf, quote=F, row.names = F, col.names = T, sep="\t")
close(gzf)

