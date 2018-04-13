require(data.table)
require(tidyverse)
DATADIR=Sys.getenv("DOX_DATA")
#DATADIR="~/gdrive/dox_data/"
genotype=fread(paste0("zcat < ",DATADIR, "genotype.txt.gz"), data.table = F, header = T)

snploc=fread(paste0("zcat < ",DATADIR, "snploc_w_rsid.txt.gz"),data.table = F) 

rownames(genotype)=genotype$snpid
genotype$snpid=NULL

genotype=as.matrix(genotype)

if (F) {
  errorCovariance=get_relatedness("../data/addSNP.coef.3671", unique(sample_anno$findiv))
  saveRDS( errorCovariance, file="../data/error_covariance.Rds" )
} else { errorCovariance = readRDS( paste0(DATADIR, "error_covariance.Rds")  ) }

input <- read.delim( paste0(DATADIR, "counts_log_cpm.txt.gz") )

anno <- read.delim( paste0(DATADIR, "sample_annotation.txt") , stringsAsFactors = F)

sample_anno=read.table( paste0(DATADIR, "annotation.txt") , header=T, stringsAsFactors = F) %>% rename(findiv=dbgap)

old_sample_anno=read.table( paste0("../data/annotation_findiv.txt") , header=T, stringsAsFactors = F)

# mapping from cell-line ID to individual
findiv=sample_anno$findiv
names(findiv)=sample_anno$cell_line
stopifnot(is.character(anno$individual))

colnames(input)=findiv[anno$individual]
findiv[ findiv=="7440_4ce2" ]="3e07_41cd" # fix sample mislabelling

anno$findiv=as.character(findiv[anno$individual])

df=data.frame( cell_line=sample_anno$cell_line, findiv=as.character(old_sample_anno$findiv), dbgap=sample_anno$findiv, stringsAsFactors = F) %>% distinct()

colnames(genotype) = data.frame(findiv=colnames(genotype), stringsAsFactors = F) %>% left_join(df, by="findiv") %>% .$dbgap

genotype=as.data.frame(genotype, stringsAsFactors = F)
genotype=cbind(snpid=rownames(genotype),genotype)

gz=gzfile(paste0(DATADIR,"genotype.txt.gz"),"w")
write.table(genotype, gz, quote=F,row.names=F, col.names=T)
close(gz)

ind_conc=fread("zcat < ../data/leafcutter_anno.txt.gz", data.table = F)
ind_conc %>% mutate(findiv=as.character(findiv)) %>% left_join(df, by="findiv") %>% select(findiv=dbgap, conc=conc) %>% write.table("../data/leafcutter_anno.txt",quote=F,row.names=F,col.names=T)
