require(data.table)
require(dplyr)
require(foreach)
source("utils.R")
require(gridExtra)
require(rstan)
# Does allelic imbalance support the response eQTLs we find? 

theme_bw(base_size = 14)

dat=fread(paste0("zcat < ",DATADIR,"/ase.txt.gz"), data.table=F) %>%
  mutate(snp=paste(chr,pos,sep=":"))  %>% 
  separate(sample, c("cell_line","cond"), sep="_")

geneloc=read.table(paste0(DATADIR,"genelocGRCh38.txt.gz"),header=T,stringsAsFactors = F)  %>% 
  filter(chr==which_chr)

bb_stan=stan_model("bb.stan")

# eqtl = read_qtls("~/gdrive/dox_data/panama_qq_boot_1e+05/")  %>%
#   mutate(p=p_interact) %>%
#   bonferroni() %>% 
#   left_join(snploc, by=c(cis_snp="snpid")) 

eqtl = 

shape_values = c(0:9, letters, LETTERS)

cell_line_to_dbgap = sample_anno %>% select(cell_line, dbgap) %>% distinct()

anno_findiv=read.table(paste0(DATADIR,"annotation_findiv.txt"), header=T, stringsAsFactors = F) %>% 
  select(cell_line, findiv) %>%
  distinct() %>% 
  left_join(sample_anno %>%
              select(cell_line, dbgap) %>%
              distinct(), 
            by="cell_line") %>%
  mutate(findiv=as.character(findiv))

phased_types=c("0|0","0|1","1|0","1|1")
phased_hets=c("0|1","1|0")


pdf("../figures/ai_support_panama.pdf",width=12,height=6)
foreach(which_chr=paste0("chr",1:22), .combine = c, .errorhandling = "stop") %do% {
  print(which_chr)
  
  top_hits = eqtl %>% filter( chr==which_chr, q < 0.1 ) 
  
  if (nrow(top_hits)==0) return(NULL)
  
  phased=fread(paste0("zcat < ", DATADIR,"phased_genotypes/dox-hg38-",which_chr,".vcf.gz"), data.table=F, skip = "#CHROM", header=T)
  #phased=fread("zcat < ../data/dox-hg38.vcf.gz", data.table=F, skip = "#CHROM", header=T)
  
  colnames(phased)[1]="CHROM"
  #phased=phased %>% filter(CHROM=="chr22")
  
  dat_chr=dat %>% filter( chr==which_chr ) %>% 
    left_join( cell_line_to_dbgap , by="cell_line" ) %>% filter( dbgap != "7440_4ce2"  ) %>%
    mutate( pos_alt=paste(pos, alt, sep="_") )
  #ase_pos=unique( dat_chr$pos )
  #length( intersect( ase_pos, phased$POS ) ) / length( ase_pos ) # 90%
  
  phased = phased %>% unite( pos_alt, POS, ALT, sep="_", remove=F) %>% 
    distinct(pos_alt, .keep_all = TRUE) 
  
  rownames( phased )=phased$pos_alt
  
  colnames(phased)[11:ncol(phased)] = 
    data.frame(findiv=colnames(phased)[11:ncol(phased)], stringsAsFactors = F) %>% 
    left_join(anno_findiv, by="findiv") %>% 
    .$dbgap
  
  dat_chr = dat_chr %>% 
    mutate( pos_alt=paste(pos, alt, sep="_")  ) %>% 
    filter( pos_alt %in% rownames(phased) ) %>% 
    mutate( geno=phased[ cbind(as.character(pos_alt), as.character( dbgap )) ] )

  # chr1 i=3 p=0.025
  foreach(i=seq_len(nrow(top_hits)), .errorhandling = "stop") %do% {
    top_hit = top_hits[i,]
    
    gene_meta = geneloc %>% filter(geneid == top_hit$gene)
    if (nrow(gene_meta)==0) return(NULL)
    
    #phased_exonic = phased %>% filter( POS > gene_meta$left, POS < gene_meta$right )
    
    ase_dat = dat_chr %>% filter( pos > gene_meta$left, pos < gene_meta$right ) 
    
    allelic_count_total=sum(ase_dat$y+ase_dat$r)
    if (allelic_count_total < 5000) {
      cat("Allelic total count only ", allelic_count_total,", skipping...\n")
      return(NULL)
    }
    
    reg_geno = phased %>% filter( POS == top_hit$pos )
    
    ase_dat$reg_geno=reg_geno[as.character(ase_dat$dbgap)] %>% as.matrix %>% as.character
    
    ge=anno %>% 
      mutate(geno=reg_geno[as.character(dbgap)] %>% 
               as.matrix %>%
               as.character, 
             y=as.numeric(input[top_hit$gene,])) %>% 
      separate(geno, c("allele1","allele2"), "[|]") %>% 
      mutate(allele1=as.numeric(allele1), 
             allele2=as.numeric(allele2), 
             geno=factor(allele1+allele2))
    
    ge_plot = ge %>% filter(!is.na(geno)) %>% ggplot(aes(as.factor(conc), 2^(y), col=geno)) + geom_boxplot() + ggtitle(paste("Gene:",top_hit$gene,"SNP:",top_hit$RSID)) + ylab("Expression (cpm)") + xlab("Dox concentration") + expand_limits(y = 0) 
    #ggsave("../figures/example_ai_vs_lrt.pdf",height=5,width=5)
    
    to_plot = ase_dat %>%
      filter( geno %in% phased_hets, reg_geno %in% phased_types ) %>% 
      # filter( geno %in% phased_hets, reg_geno %in% phased_hets ) %>% 
      mutate( in_phase=geno == reg_geno ) %>% 
      filter( (r+y) > 0 ) %>% 
      mutate( ind=as.factor(dbgap), snp=as.factor(pos), coverage=r+y, ar = r/coverage, car=ifelse(in_phase,ar,1-ar) )
    
    if (nrow(to_plot)==0) return(NULL)
  
    pv=anova( lm(car ~ 1, data=to_plot, weights = sqrt(r+y)), lm(car ~ cond, data=to_plot, weights = sqrt(r+y) ) )[2,6]
    
    levels(to_plot$snp)
    #to_plot %>% ggplot(aes( cond,car,col=ind,size=coverage,shape=snp)) + geom_point(position = position_jitter(width = 0.3, height = 0)) + ylim(0,1) + ggtitle(paste0(top_hit$gene," p=",format(pv,digits=3))) + scale_shape_manual(values=shape_values[seq_along(levels(to_plot$snp))] )  
   
   x_df = to_plot %>% mutate( het_x=ifelse(reg_geno %in% phased_hets, ifelse(in_phase,1,-1), 0) , cond=factor(cond))
   x_full = model.matrix( ~ het_x + cond:het_x, data=x_df )
   stan_dat = list(N=nrow(x_full), P=ncol(x_full), x=x_full, ys=to_plot$y, ns=to_plot$y + to_plot$r, concShape=1.001, concRate=0.001)
   fit_full=optimizing(bb_stan, data=stan_dat)
   
   x_null = model.matrix( ~ het_x, data=x_df )
   stan_dat = list(N=nrow(x_null), P=ncol(x_null), x=x_null, ys=to_plot$y, ns=to_plot$y + to_plot$r, concShape=1.001, concRate=0.001)
   fit_null=optimizing(bb_stan, data=stan_dat)
   
   df=ncol(x_full) - ncol(x_null)
   lrt=2.0*(fit_full$value - fit_null$value)
   #-pchisq(lrt, df, lower.tail = F, log.p = T)/log(10)
   pv=pchisq(lrt, df, lower.tail = F)
   
   ase_plot = 
     to_plot %>% mutate(het=reg_geno %in% phased_hets) %>% ggplot(aes( cond,car,size=coverage,shape=snp,col=het)) + geom_point() + ylim(0,1) + ggtitle(paste0(top_hit$gene," p=",format(pv,digits=3))) + scale_shape_manual(values=shape_values[seq_along(levels(to_plot$snp))], guide=F )   + geom_line(aes(group=interaction(ind,snp)),size=1,alpha=0.5)  + xlab("Dox concentration") + ylab("Phased allelic ratio")
   
   print("Plotting")
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
