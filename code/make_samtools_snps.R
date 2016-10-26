# Make a common SNP table which samtools mpileup will use. 
require(data.table)

# Downloaded from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/snp147Common.txt.gz
snps=fread("zcat < ~/scailscratch/Homo_sapiens/NCBI/GRCh38/snp147Common.txt.gz")

setDF(snps)
snps_small=snps[,2:3]

colnames(snps_small)=c("chr","pos")

# Not sure what these SNPs with position=0 (only 3 of them) represent, but samtools mpileup
# doesn't like them. 
snps_small=snps_small[ snps_small$pos > 0, ]

snps_small$pos = snps_small$pos + 1

write.table(snps_small, file="~/scailscratch/Homo_sapiens/NCBI/GRCh38/snp147.txt", col.names = F, row.names = F, quote=F)

# This was code to check for a reasonable overlap between hets from RNA-seq vs common snps. 
if (FALSE) {

  snps=read.table("~/scailscratch/Homo_sapiens/NCBI/GRCh38/snp147.txt", header=F, sep=" ", stringsAsFactors=F)
  
  test=fread("zcat < ~/scailscratch/dox/temp")
  
  setDF(test)
  colnames(test)=c("chr","pos","ref","alt","ref_allele","alt_allele")
  
  colnames(snps)=c("chr","pos")
  
  chr="chr1"
  test_chr=test[test$chr==chr,]
  snps_chr=snps[snps$chr==chr,]
  
  inte=intersect(test_chr$pos, snps_chr$pos)
  
}
  
