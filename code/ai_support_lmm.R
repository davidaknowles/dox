require(data.table)
require(dplyr)
require(foreach)
source("utils.R")

theme_bw(base_size = 14)

dat=fread("zcat < ../data/ase.txt.gz", data.table=F) %>%
  mutate(snp=paste(chr,pos,sep=":"))

cl_cond=do.call( rbind, strsplit( dat$sample, "_" ) )
dat$cell_line=cl_cond[,1]
dat$cond=cl_cond[,2]

sample_anno=read.table("../data/annotation.txt", header=T, stringsAsFactors = F) %>%
  select( cell_line, findiv ) %>%
  distinct()

geneloc=read.table("../data/genelocGRCh38.txt.gz",header=T,stringsAsFactors = F)
rownames(geneloc)=geneloc$geneid

input <- read.delim("../data/counts_log_cpm.txt.gz", check.names = F)
# corresponds to cols of input
anno <- read.delim("../data/sample_annotation.txt", stringsAsFactors = F)  %>% 
  left_join(sample_anno, by=c(individual="cell_line")) %>%
  mutate(findiv=as.character(findiv))

snploc=read.table("../data/snploc.txt.gz",header=T,stringsAsFactors = F)
#snploc$chrpos=with(snploc, paste(chr,pos,sep="_"))

lrt_res=read.table("../data/lrt-summary.txt.gz", stringsAsFactors = F, header=T) %>% 
  left_join(snploc, by=c("snp"="snpid") ) %>% 
  mutate( q=p.adjust( pmin( p, 1), method = "BH" ) )

shape_values = c(0:9, letters, LETTERS)

pdf("../figures/ai_support_lrt.pdf",width=12,height=6)
foreach(which_chr=paste0("chr",1:22), .combine = c, .errorhandling = "stop") %do% {
  print(which_chr)
  
  top_hits = lrt_res %>% filter( chr==which_chr, q < 0.1 ) 
  
  if (nrow(top_hits)==0) return(NULL)
  
  phased=fread(paste0("zcat < ../data/phased_genotypes/dox-hg38-",which_chr,".vcf.gz"), data.table=F, skip = "#CHROM", header=T)
  #phased=fread("zcat < ../data/dox-hg38.vcf.gz", data.table=F, skip = "#CHROM", header=T)
  
  colnames(phased)[1]="CHROM"
  #phased=phased %>% filter(CHROM=="chr22")
  
  dat_chr=dat %>% filter( chr==which_chr ) %>% 
    left_join( sample_anno , by="cell_line" ) %>% filter( findiv != 160001  ) %>%
    mutate( pos_alt=paste(pos, alt, sep="_") )
  #ase_pos=unique( dat_chr$pos )
  #length( intersect( ase_pos, phased$POS ) ) / length( ase_pos ) # 90%
  
  phased = phased %>% mutate( pos_alt=paste(POS, ALT, sep="_") ) %>% 
    distinct(pos_alt, .keep_all = TRUE) 
  
  rownames( phased )=phased$pos_alt
  
  dat_chr = dat_chr %>% 
    mutate( pos_alt=paste(pos, alt, sep="_")  ) %>% 
    filter( pos_alt %in% rownames(phased) ) %>% 
    mutate( geno=phased[ cbind(as.character(pos_alt), as.character( findiv )) ] )

  # chr1 i=3 p=0.025
  foreach(i=seq_len(nrow(top_hits)), .errorhandling = "stop") %do% {
    top_hit = top_hits[i,]
    
    gene_meta = geneloc %>% filter(geneid == top_hit$gene)
    
    #phased_exonic = phased %>% filter( POS > gene_meta$left, POS < gene_meta$right )
    
    ase_dat = dat_chr %>% filter( pos > gene_meta$left, pos < gene_meta$right ) 
    
    if (sum(ase_dat$y+ase_dat$r) < 5000) return(NULL)
    
    reg_geno = phased %>% filter( POS == top_hit$pos )
  
    ase_dat$reg_geno=reg_geno[as.character(ase_dat$findiv)] %>% as.matrix %>% as.character
    
    ge=anno %>% mutate(geno=reg_geno[as.character(findiv)] %>% as.matrix %>% as.character, y=as.numeric(input[top_hit$gene,]))
    ge$geno = foreach(s=strsplit(ge$geno,"|",fixed=T)) %do% { (2-as.numeric(s)) %>% sum } %>% unlist %>% factor
    
    ge_plot = ge %>% filter(!is.na(geno)) %>% ggplot(aes(as.factor(conc), 2^(y), col=geno)) + geom_boxplot() + ggtitle(paste("Gene:",top_hit$gene,"SNP:",top_hit$snp)) + ylab("Expression (cpm)") + xlab("Dox concentration") + expand_limits(y = 0) 
    ggsave("../figures/example_ai_vs_lrt.pdf",height=5,width=5)
    
    phased_types=c("0|0","0|1","1|0","1|1")
    to_plot = ase_dat %>%
      # filter( geno %in% phased_types, reg_geno %in% phased_types ) %>% 
      filter( geno %in% c("0|1","1|0"), reg_geno %in% c("0|1","1|0") ) %>% 
      mutate( in_phase=geno == reg_geno ) %>% 
      filter( (r+y) > 0 ) %>% 
      mutate( ind=as.factor(findiv), snp=as.factor(pos), coverage=r+y, ar = r/coverage, car=ifelse(in_phase,ar,1-ar) )
    
    if (nrow(to_plot)==0) return(NULL)
  
    pv=anova( lm(car ~ 1, data=to_plot, weights = sqrt(r+y)), lm(car ~ cond, data=to_plot, weights = sqrt(r+y) ) )[2,6]
  
    levels(to_plot$snp)
    #to_plot %>% ggplot(aes( cond,car,col=ind,size=coverage,shape=snp)) + geom_point(position = position_jitter(width = 0.3, height = 0)) + ylim(0,1) + ggtitle(paste0(top_hit$gene," p=",format(pv,digits=3))) + scale_shape_manual(values=shape_values[seq_along(levels(to_plot$snp))] )  
    
   ase_plot = to_plot %>% ggplot(aes( cond,car,col=ind,size=coverage,shape=snp)) + geom_point() + ylim(0,1) + ggtitle(paste0(top_hit$gene," p=",format(pv,digits=3))) + scale_shape_manual(values=shape_values[seq_along(levels(to_plot$snp))] )   + geom_line(aes(group=interaction(ind,snp) ),size=0.5,alpha=0.5)  + xlab("Dox concentration") + ylab("Phased allelic ratio")
   
   print( grid.arrange(ge_plot, ase_plot, nrow=1 ) )
   # %>% print
    # geom_dotplot(binaxis = "y", stackdir = "center", binwidth=0.01, alpha=.5)
   NULL
  }
  # ase_dat = ase_dat %>% mutate( load=as.numeric(substr(geno,1,1)) + as.numeric(substr(geno,3,3)) )
  
  # ase_dat %>% filter( (r+y) > 0 ) %>%  mutate( ar = r/(r+y) ) %>% ggplot(aes(cond,ar,fill=as.factor(load))) + geom_violin() + geom_boxplot(outlier.shape = NA, width=0.4) + geom_dotplot(binaxis = "y", stackdir = "center", binwidth=0.01, alpha=.5)
  
  # ase_dat %>% filter( (r+y) > 20 ) %>%  mutate( ar = r/(r+y) ) %>% ggplot(aes(cond,ar,fill=as.factor(load))) + geom_boxplot()

}
dev.off()
