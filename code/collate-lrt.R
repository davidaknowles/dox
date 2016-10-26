require(ggplot2)
require(doMC)

DATADIR="~/scailscratch/dox/"
library("dplyr")
library(data.table)
source("utils.R")

genotype=fread(paste0("zcat < ",DATADIR,"genotype.txt.gz"))
setDF(genotype)
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

pdf("../figures/lrt-iqtl.pdf",width=7,height=5)
res= foreach(chrom=chroms, .combine = c) %do% {
  print(chrom)
  fn=paste0(DATADIR,"lrt_",chrom,".RData")
  if (!file.exists(fn)) {
      cat("No file:",fn,"\n")
      return(NA)
  }
  load(fn)
  foreach(gene=names(results), .combine = c) %do% {
      result=results[[gene]]
      ps=foreach(cis_snp=names(result), .combine=c) %do% {
          res=result[[cis_snp]]
          if (!("error" %in% class(res))) res$p else NA
      }

    bfp=min(ps,na.rm=T) * sum(!is.na(ps))
    if (!is.na(bfp)) if (bfp < 0.0005) {
      print(gene)
      cis_snp=names(result)[which.min(ps)]
      y=input[gene,]
      geno=genotype[cis_snp,anno$findiv]
      print(ggplot(data.frame(y=y, geno=as.factor(geno), conc=anno$conc), aes(as.factor(conc), y, col=geno)) + geom_boxplot() + ggtitle(paste("Gene:",gene,"SNP:",cis_snp)) + ylab("Expression") + xlab("Dox concentration") + theme_bw(base_size=16))
      print(c(cis_snp, gene))
    }

    bfp
  } 
} 
dev.off() 

q=p.adjust(res, method="BH")
sum(q<.05,na.rm = T)
# q<0.05 ~ p<0.0005
