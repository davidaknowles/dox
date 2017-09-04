DATADIR=Sys.getenv("DOX_DATA")
#DATADIR="~/gdrive/dox/"
genotype=fread(paste0("zcat < ",DATADIR, "genotype.txt.gz"), data.table = F, header = T)

snploc=fread(paste0("zcat < ",DATADIR, "snploc.txt.gz"),data.table = F) 

rownames(genotype)=genotype$snpid
genotype$snpid=NULL

genotype=as.matrix(genotype)

if (F) {
  errorCovariance=get_relatedness("../data/addSNP.coef.3671", unique(sample_anno$findiv))
  saveRDS( errorCovariance, file="../data/error_covariance.Rds" )
} else { errorCovariance = readRDS( paste0("zcat < ",DATADIR, "error_covariance.Rds")  ) }


input <- read.delim( paste0("zcat < ",DATADIR, "counts_log_cpm.txt.gz") )

anno <- read.delim( paste0("zcat < ",DATADIR, "sample_annotation.txt") , stringsAsFactors = F)

sample_anno=read.table( paste0("zcat < ",DATADIR, "annotation.txt") , header=T, stringsAsFactors = F)

# mapping from cell-line ID to individual
findiv=sample_anno$findiv
names(findiv)=sample_anno$cell_line
stopifnot(is.character(anno$individual))

colnames(input)=findiv[anno$individual]
findiv[ findiv==160001 ]=106411

anno$findiv=as.character(findiv[anno$individual])