require(data.table)

DATADIR=Sys.getenv("DOX_DATA")
#DATADIR="~/gdrive/dox_data/"
genotype=fread(paste0("zcat < ",DATADIR, "genotype.txt.gz"), data.table = F, header = T)
rownames(genotype)=genotype$snpid
genotype$snpid=NULL
genotype=as.matrix(genotype)

anno_findiv=read.table( "../data/annotation_findiv.txt" , header=T, stringsAsFactors = F)  %>% 
  select(cell_line, findiv) %>% 
  distinct()

sample_anno=read.table( paste0(DATADIR, "annotation.txt") , header=T, stringsAsFactors = F)

anno_findiv = sample_anno %>% select(cell_line, dbgap) %>% distinct() %>% left_join(anno_findiv, by="cell_line")

colnames(genotype)=data.frame(findiv=as.integer(colnames(genotype))) %>% left_join(anno_findiv, by="findiv") %>% .$dbgap

gz=gzfile(paste0(DATADIR, "genotype_dbgap.txt.gz"),"w")
genotype %>% write.table(gz, quote=F, row.names=T, col.names=T, sep="\t")
close(gz)
