source("utils.R")
require(doMC)
require(MatrixEQTL)
require(irlba)
require(Matrix)

args=if (interactive()) c("0.0","10") else commandArgs(trailingOnly = T)

use_pedigree=T

DATADIR="~/scailscratch/dox/"

#covariate_file=paste0(DATADIR,"covariates.txt")
covariate_file=paste0(DATADIR,"covariates_with_PC.txt")

# Load gene expression data
input <- read.delim("../data/counts_log_cpm.txt.gz")
anno <- read.delim("../data/sample_annotation.txt", stringsAsFactors = F)

conc_to_run=as.numeric(args[1])
num_PCs_to_remove=as.integer(args[2])

sample_anno=read.table("../data/annotation.txt", header=T, stringsAsFactors = F)

# mapping from cell-line ID to individual
findiv=sample_anno$findiv
names(findiv)=sample_anno$cell_line
stopifnot(is.character(anno$individual))

colnames(input)=findiv[anno$individual]

input=remove_PCs(input, num_PCs_to_remove)

input=quantile_normalize(input)

# Write GE to file
ge_conc=qnorm_input[,anno$conc==conc_to_run]
ge_filename=paste0(DATADIR, "ge", conc_to_run, "_PC", num_PCs_to_remove, ".txt")
write.table(ge_conc, file=ge_filename, sep="\t", quote=F)


# Load IBD based relatedness matrix
if (use_pedigree) {
  errorCovariance=get_relatedness("../data/addSNP.coef.3671", unique(sample_anno$findiv))
}

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = " ";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(paste0(DATADIR, "genotype.txt"));

## Look at the distribution of MAF
# hist(maf[maf<0.1],seq(0,0.1,length.out=100))

cat('SNPs before filtering:',nrow(snps))
# snps$RowReorderSimple(maf>0.05);
#snps$RowReorder(maf>0.05);
cat('SNPs before filtering:',nrow(snps))

## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(ge_filename)

# Line up GE and genotype samples
shared_ind=intersect(gene$columnNames, snps$columnNames)

temp=1:length(gene$columnNames)
names(temp)=gene$columnNames
gene$ColumnSubsample(temp[shared_ind])

temp=1:length(snps$columnNames)
names(temp)=snps$columnNames
snps$ColumnSubsample(temp[shared_ind])

# Check MAF _after_ subsampling individuals
maf=foreach(sl=seq_len(length(snps)), .combine=c) %do% {
  m = rowMeans(snps[[sl]],na.rm=TRUE)/2;
  pmin(m,1-m)
}
length(maf)
snps$RowReorderSimple(maf>0.05)

if (use_pedigree) {
errorCovariance=errorCovariance[shared_ind,shared_ind]
} else errorCovariance=numeric()


cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$LoadFile(covariate_file);

temp=1:length(cvrt$columnNames)
names(temp)=cvrt$columnNames
cvrt$ColumnSubsample(temp[shared_ind])

## Run the analysis
snpspos = read.table(paste0(DATADIR, "snploc.txt"), header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(paste0(DATADIR, "genelocGRCh38.txt"), header = TRUE, stringsAsFactors = FALSE);

suffix=paste0(conc_to_run, "_PC", num_PCs_to_remove)

me = Matrix_eQTL_main(
		snps = snps, 
		gene = gene, 
		cvrt = cvrt,
		output_file_name     = NULL,
		pvOutputThreshold     = 0,
		useModel = modelLINEAR, 
		errorCovariance = errorCovariance, 
		verbose = TRUE, 
		output_file_name.cis = paste0(DATADIR, "results_pc_", suffix, ".txt"),
		pvOutputThreshold.cis = 1,
		snpspos = snpspos, 
		genepos = genepos,
		cisDist = cisDist,
		pvalue.hist = TRUE,
		min.pv.by.genesnp = TRUE,
		noFDRsaveMemory = FALSE);

## Results:
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');

qtls=me$cis$eqtls
colnames(qtls)[1]="SNP"
colnames(qtls)[4]="p.value"
require(dplyr)
myq=qtls %>% group_by(gene) %>% summarize( SNP=SNP[which.min(p.value)], p=min(p.value)*length(p.value))
myq$p[myq$p>1]=1
myq$q=p.adjust(myq$p,method="BH")
gzf=gzfile(paste0(basedir,"summary",suffix),"w")
write.table(as.data.frame(myq), gzf, quote=F, sep="\t", row.names=F)
close(gzf)
nsig=sum(myq$q<0.05,na.rm=T)
print(nsig)
nsig
