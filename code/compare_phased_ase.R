require(data.table)
require(dplyr)
require(foreach)

dat=fread("zcat < ../data/ase_full.txt.gz", data.table=F)

dat$snp=with(dat, paste(chr,pos,sep=":"))

which_chr="chr1"

phased=fread(paste0("zcat < ../data/phased_genotypes/dox-hg38-",which_chr,".vcf.gz"), data.table=F, skip = "#CHROM", header=T)
#phased=fread("zcat < ../data/dox-hg38.vcf.gz", data.table=F, skip = "#CHROM", header=T)

colnames(phased)[1]="CHROM"
#phased=phased %>% filter(CHROM=="chr22")

dat_chr=dat %>% filter( chr==which_chr )
samp_cond=do.call( rbind, strsplit( dat_chr$sample, "_" ) )
dat_chr$cell_line=samp_cond[,1]
dat_chr$cond=samp_cond[,2]
dat_chr$sample=NULL

ase_pos=unique( dat_chr$pos )
length( intersect( ase_pos, phased$POS ) ) / length( ase_pos ) # 90%

sample_anno=read.table("../data/annotation.txt", header=T, stringsAsFactors = F)

sample_anno=sample_anno  %>% select( cell_line, dbgap ) %>% distinct()

dat_chr = dat_chr %>% left_join( sample_anno , by="cell_line" )

phased = phased %>% mutate( pos_alt=paste(POS, ALT, sep="_") ) %>% distinct(pos_alt, .keep_all = TRUE)

rownames( phased )=phased$pos_alt

#verify=read.table("../data/verify.txt", stringsAsFactors = F, header=T)
#colnames(verify)[1:2]=c("best","dbgap")

#dat_chr = dat_chr %>% left_join(verify %>% select( dbgap, best ), by="dbgap")
#dat_chr$dbgap=dat_chr$best

dat_chr = dat_chr %>% mutate( pos_alt=paste(pos, alt, sep="_")  ) %>% filter( pos_alt %in% rownames(phased) )

dat_chr=dat_chr %>% mutate ( geno=phased[ cbind(as.character(pos_alt), as.character( dbgap )) ] ) %>% mutate( load=as.numeric(substr(geno,1,1)) + as.numeric(substr(geno,3,3)) )

pdf("../figures/7440_4ce2_is_bad.pdf",height=6,width=8)
foreach(which_ind=unique(dat_chr$dbgap)) %do% {
  ggplot( dat_chr %>% filter( (r+y)>10 , dbgap==which_ind, !is.na(load)) , aes(y/(r+y), fill=as.factor(load) )) + facet_wrap( ~sample )  + geom_histogram(alpha=0.5, position="identity")  + scale_fill_discrete("Genotype") + xlab("Allelic ratio") + ggtitle(which_ind)
  }
dev.off()

dat_chr %>% 
  filter( (r+y)>10 , dbgap=="7440_4ce2" ) %>% 
  mutate ( geno=phased[ cbind(as.character(pos_alt), "3e07_41cd") ] ) %>% 
  mutate( load=as.numeric(substr(geno,1,1)) + as.numeric(substr(geno,3,3)) ) %>% 
  filter( !is.na(load)) %>% 
  ggplot( aes(y/(r+y), fill=as.factor(load) )) + facet_wrap( ~sample )  + geom_histogram(alpha=0.5, position="identity")  + scale_fill_discrete("Genotype") + xlab("Allelic ratio") + ggtitle("7440_4ce2 RNA and 3e07_41cd genotype")

add_counts = dat_chr %>% group_by( pos_alt, ref, alt, dbgap ) %>% summarize( r=sum(r), y=sum(y) )

theme_set(theme_bw(base_size = 14))
require(tidyr)
require(magrittr)
# 3e07_41cd
comp_dat= dat_chr %>% filter( (r+y)>30 , r>2, y>2 ) %>% mutate( ar=y/(r+y) ) %>% filter( abs(ar-0.5)<0.4 ) %>% select(dbgap, ar, pos_alt, cond) %>% spread( dbgap, ar )
colnames(comp_dat)=make.names(colnames(comp_dat))
comp_dat %>% ggplot( aes(x=X3e07_41cd,y=X7440_4ce2) ) + geom_point(alpha=.3) + xlim(0,1) + ylim(0,1)  + facet_wrap( ~cond )

comp_dat %>% group_by(cond) %>% summarise(r2=cor(X111011,X7440_4ce2,use="pairwise")^2,p=cor.test(X111011,X7440_4ce2,use="pairwise")$p.value)
comp_dat %>% group_by(cond) %>% summarise(r2=cor(X3e07_41cd,X7440_4ce2,use="pairwise")^2,p=cor.test(X3e07_41cd,X7440_4ce2,use="pairwise")$p.value)


counts <- read.delim("../data/counts_log_cpm.txt.gz") 
#counts=read.table("../data/counts.txt.gz", check.names = F)
anno <- read.delim("../data/sample_annotation.txt", stringsAsFactors = F)
anno = anno %>% left_join( sample_anno , by=c("individual"="cell_line" ) )


inds=unique(anno$dbgap)
count_cors=foreach(co=unique(anno$conc)) %dopar% {
  foreach(ind1=inds, .combine = rbind) %do% {
    foreach(ind2=inds, .combine = c) %do% {
      replicates=counts[,which( anno$conc==co & anno$dbgap %in% c(ind1,ind2) ), drop=F ]
      if (ncol(replicates)<2) return(NA)
      cor(replicates[,1], replicates[,2], use="pairwise")
    }
  }
}

foreach(i=1:5) %do%
{ dimnames(count_cors[[i]])=list(inds,inds) }

count_cors[[1]][ lower.tri(count_cors[[1]]) ]=NA
count_cors[[1]][ "7440_4ce2", "3e07_41cd" ]

foreach(i=1:5) %do%
{ dimnames(count_cors[[i]])=list(inds,inds)
 mean(as.numeric(count_cors[[i]]) > count_cors[[i]][ "7440_4ce2", "3e07_41cd" ], na.rm=T) }

hist(as.numeric(count_cors[[1]]))

ind1="7440_4ce2"
ind2="3e07_41cd"
require(gridExtra)
do.call(grid.arrange,c(
foreach(co=unique(anno$conc)) %do% {
  replicates=counts[,which( anno$conc==co & anno$dbgap %in% c(ind1,ind2) ) ]
  colnames(replicates)=c("ind1","ind2")
  ggplot( replicates, aes(ind1, ind2)) + geom_point() + scale_x_log10(limits=c(1,1e6)) + scale_y_log10(limits=c(1,1e6)) + ggtitle(paste("Dox conc",co)) + geom_abline(intercept = 0, slope=1, col="red")
}, nrow=2 ))
#phased_mat=phased[,10:ncol(phased)]
#colnames(phased_mat)=sample(colnames(phased_mat))

#add_counts$geno=phased_mat[ cbind(as.character(add_counts$pos_alt), as.character( add_counts$dbgap )) ]

add_counts$REF=phased[ as.character(add_counts$pos_alt), "REF" ]
add_counts$ALT=phased[ as.character(add_counts$pos_alt), "ALT" ]

with(add_counts, mean(ref==REF)) # 99.8%
with(add_counts, mean(alt==ALT)) # 99.4% # now only 56

table( add_counts$geno )

add_counts = add_counts[ grep("\\.",add_counts$geno,invert = T), ]
add_counts$hap1=as.numeric(substr(add_counts$geno,1,1))
add_counts$hap2=as.numeric(substr(add_counts$geno,3,3))
add_counts$load=with(add_counts, hap1+hap2)
add_counts$het=with(add_counts, (hap1+hap2)==1)
add_counts$het=with(add_counts, (hap1+hap2)==1)

ggplot( add_counts %>% filter( (r+y)>10 , dbgap=="7440_4ce2") , aes(y/(r+y), fill=as.factor(load) ))  + geom_histogram(alpha=0.5, position="identity") + theme(legend.position = c(.8,.8)) + scale_fill_discrete("Genotype") + xlab("Allelic ratio")


problem=dat_chr %>% filter( (r+y)>10 , dbgap=="7440_4ce2", !is.na(load) ) %>% mutate( ar=y/(r+y) ) 
samples=unique(problem$sample)

d=dcast(problem, pos_alt ~ sample, value.var="ar")

for (i in seq_len(length(samples)-1))
  for (j in (i+1):length(samples))
  {
    plot(d[,1+i],d[,1+j])
  }

