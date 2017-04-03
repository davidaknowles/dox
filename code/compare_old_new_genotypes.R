require(data.table)
which_chr="chr22"
phased=fread(paste0("zcat < ../data/phased_genotypes/dox-hg38-",which_chr,".vcf.gz"), data.table=F, skip = "#CHROM", header=T)
genotype=fread("zcat < ../data/genotype.txt.gz", data.table = F, header = T)
rownames(genotype)=genotype$snpid
genotype$snpid=NULL

snploc=read.table("../data/snploc.txt.gz",header=T,stringsAsFactors = F)

genotype=genotype[snploc$chr==which_chr,] # 1% missing, 30176 x 46
snploc_chr=snploc[snploc$chr==which_chr,]

phased=phased[!duplicated(phased$POS),]

require(foreach)
require(doMC)
registerDoMC(7)
phased_load=foreach(i=10:ncol(phased), .combine = cbind) %dopar% {
  foreach(j=c(1,3), .combine = "+") %do% { as.numeric( substr( phased[,i], j, j ) ) }
}

rownames(phased_load)=phased$POS

temp=phased_load[ as.character(snploc_chr$pos), ]

mean(is.na(genotype))
mean(is.na(temp))

missing_per_ind=colMeans( is.na( phased_load ) )
hist(missing_per_ind)

missing_per_snp=rowSums( is.na( phased_load ) )
hist(missing_per_snp)

sum(missing_per_snp <= 2) # only 39871 compared to 30176 genotyped, so probably not worth running iQTLs on imputed SNPs? 

to_keep=missing_per_ind < 0.8
phased_load=phased_load[,to_keep]

to_keep= rowSums(is.na(phased_load)) <= 10
phased_load=phased_load[to_keep,]


