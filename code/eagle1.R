require(data.table)
require(dplyr)
require(foreach)
source("utils.R")
source("load_data.R")
require(rstan)
require(doMC)
registerDoMC(16)

cisdist=1e5

if (interactive()) {
  which_chr="chr22"
} else {
  which_chr=commandArgs(trailingOnly = T)[1]
}

cat("Chrom:",which_chr,"\n")

dat=fread("zcat < ../data/ase.txt.gz", data.table=F) %>%
  mutate(snp=paste(chr,pos,sep=":"))  %>% 
  separate(sample, c("cell_line","cond"), sep="_")

geneloc=read.table(paste0(DATADIR,"genelocGRCh38.txt.gz"),header=T,stringsAsFactors = F)
rownames(geneloc)=geneloc$geneid

bb_stan=stan_model("bb.stan")

anno_findiv=read.table("../data/annotation_findiv.txt", header=T, stringsAsFactors = F) %>% 
  select(cell_line, findiv) %>%
  distinct() %>% 
  left_join(sample_anno %>%
              select(cell_line, dbgap) %>%
              distinct(), 
            by="cell_line") %>%
  mutate(findiv=as.character(findiv))

phased_types=c("0|0","0|1","1|0","1|1")
phased_hets=c("0|1","1|0")

phased=fread(paste0("zcat < ../data/phased_genotypes/dox-hg38-",which_chr,".vcf.gz"), data.table=F, skip = "#CHROM", header=T)
#phased=fread("zcat < ../data/dox-hg38.vcf.gz", data.table=F, skip = "#CHROM", header=T)

colnames(phased)[1]="CHROM"
#phased=phased %>% filter(CHROM=="chr22")

dat_chr=dat %>% filter( chr==which_chr ) %>% 
  left_join( sample_anno , by="cell_line" ) %>% filter( dbgap != "7440_4ce2"  ) %>%
  mutate( pos_alt=paste(pos, alt, sep="_") )

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

outputdir=paste0(DATADIR,"/eagle1/")
dir.create( outputdir, showWarnings = F )
checkpoint_dir=paste0(outputdir,which_chr,"/")
dir.create(checkpoint_dir, recursive = T, showWarnings = F)

foreach(gene=unique(geneloc$geneid), .errorhandling = "remove", .combine = bind_rows) %do% {

  if (!is.null(checkpoint_dir)){
    check_fn=paste0(checkpoint_dir,gene,".txt.gz")
    if (file.exists(check_fn)) {
      return(read.table(check_fn, header=T, stringsAsFactors = F, sep="\t"))
    }
  }
  gene_meta = geneloc %>% filter(geneid == gene)
  
  #phased_exonic = phased %>% filter( POS > gene_meta$left, POS < gene_meta$right )
  
  ase_dat = dat_chr %>% filter( pos > gene_meta$left, pos < gene_meta$right ) 
  
  allelic_count_total=sum(ase_dat$y+ase_dat$r)
  if (allelic_count_total < 5000) {
    cat("Allelic total count only ", allelic_count_total,", skipping...\n")
    return(NULL)
  }
  
  if (nrow(ase_dat) < 20) return(NULL)
  
  cis_snps=phased %>% 
    filter( (geneloc[gene,"left"]-cisdist) < POS, 
            (geneloc[gene,"right"]+cisdist) > POS ) %>% .$POS
  
  gene_results = foreach(snp_pos=cis_snps, .errorhandling = "remove",  .combine = bind_rows) %dopar% {
    print(snp_pos)
    reg_geno = (phased %>% filter( POS == snp_pos ))[,11:ncol(phased)] %>% as.matrix()
    
     if (nrow(reg_geno) != 1) {
       print("Skipping >biallelic site")
       return(NULL)
     }
    
    reg_geno=data.frame(dbgap=colnames(reg_geno), reg_geno=as.character(reg_geno), stringsAsFactors = F)
    
    ase_temp = ase_dat %>% inner_join(reg_geno, by="dbgap") %>%
      filter( geno %in% phased_hets, reg_geno %in% phased_types, (r+y) > 0 ) %>% 
      mutate( het_x=ifelse(reg_geno %in% phased_hets, ifelse(geno == reg_geno,1,-1), 0) , 
              cond=factor(cond))
    
    if (nrow(ase_temp) < 10) return(NULL)
    if (length(unique(ase_temp$cond)) <= 1) return(NULL) # only data for one 
    if (sum(ase_temp$het_x != 0) < 10) return(NULL) # no heterozygous regulatory SNPs
    
    x_full = if (length(unique(ase_temp$het_x)) > 1) model.matrix( ~ het_x + cond:het_x, data=ase_temp ) else model.matrix( ~ cond, data=ase_temp ) # testing the exonic SNP itself
    
     stan_dat = list(N=nrow(x_full), P=ncol(x_full), x=x_full, ys=ase_temp$y, ns=ase_temp$y + ase_temp$r, concShape=1.001, concRate=0.001)
     fit_full=if (det(t(x_full) %*% x_full) > 0) { optimizing(bb_stan, data=stan_dat)$value } else NA

     x_null = if (length(unique(ase_temp$het_x)) > 1) model.matrix( ~ het_x , data=ase_temp ) else model.matrix( ~ 1, data=ase_temp )
     stan_dat$x=x_null
     stan_dat$N=nrow(x_null)
     stan_dat$P=ncol(x_null)
     fit_null=optimizing(bb_stan, data=stan_dat)$value
     
     x_0 = model.matrix( ~ 1, data=ase_temp )
     stan_dat$x=x_0
     stan_dat$N=nrow(x_0)
     stan_dat$P=ncol(x_0)
     fit_0=optimizing(bb_stan, data=stan_dat)$value
     
     #ase_temp %>% mutate(het=reg_geno %in% phased_hets, coverage=r+y, ar = r/coverage,in_phase=geno == reg_geno,  car=ifelse(in_phase,ar,1-ar)) %>% ggplot(aes( cond,car,size=coverage,shape=snp,col=het)) + geom_point() + ylim(0,1) + ggtitle(paste0(top_hit$gene," p=",format(pv,digits=3))) + scale_shape_manual(values=shape_values[seq_along(levels(to_plot$snp))], guide=F )   + geom_line(aes(group=interaction(dbgap,snp)),size=1,alpha=0.5)  + xlab("Dox concentration") + ylab("Phased allelic ratio")
     
     df=ncol(x_full) - ncol(x_null)
     #lrt=2.0*(fit_full$value - fit_null$value)
     # nlog10p=-pchisq(lrt, df, lower.tail = F, log.p = T)/log(10)
     data.frame(gene=gene, snp=snp_pos, df=df, l0=fit_0, l1=fit_null, l2=fit_full)
  }
  
  if (!is.null(checkpoint_dir)){
    print("Saving results")
    checkpoint_file= gzfile( check_fn,"w")
    gene_results %>% format(digits=5) %>% write.table(checkpoint_file, quote = F, row.names = F, col.names = T, sep="\t")
    close(checkpoint_file)
  }
  gene_results
  
}
