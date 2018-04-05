require(data.table)

DATADIR=Sys.getenv("DOX_DATA")
#DATADIR="~/gdrive/dox_data/"
genotype=read.table(paste0(DATADIR, "genotype_dbgap.txt.gz"), stringsAsFactors = F, sep="\t", header = T, check.names = F)

snploc=fread(paste0("zcat < ",DATADIR, "snploc_w_rsid.txt.gz"),data.table = F) 

#rownames(genotype)=genotype$snpid
#genotype$snpid=NULL

genotype=as.matrix(genotype)

if (F) {
  errorCovariance=get_relatedness("../data/addSNP.coef.3671", unique(sample_anno$dbgap))
  saveRDS( errorCovariance, file="../data/error_covariance.Rds" )
} else { errorCovariance = readRDS( paste0(DATADIR, "error_covariance.Rds")  ) }

input <- read.delim( paste0(DATADIR, "counts_log_cpm.txt.gz") )

anno <- read.delim( paste0(DATADIR, "sample_annotation.txt") , stringsAsFactors = F)

sample_anno=read.table( paste0(DATADIR, "annotation.txt") , header=T, stringsAsFactors = F)

# mapping from cell-line ID to individual
dbgap=sample_anno$dbgap
names(dbgap)=sample_anno$cell_line
stopifnot(is.character(anno$individual))

colnames(input)=dbgap[anno$individual]
dbgap[ dbgap=="7440_4ce2" ]="3e07_41cd"

anno$dbgap=as.character(dbgap[anno$individual])