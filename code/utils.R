require(Matrix)

read_qtls = function(res_dir, df=4) {

  sqtl=foreach(fn=list.files(res_dir,glob2rx("chr*.txt.gz")), .combine = bind_rows) %do% {
    cat(".")
    if (fn=="chrX.txt.gz") return(NULL)
    fread(paste0("zcat < ",res_dir,fn), data.table=F)
  }
  
  df=4
  sqtl %>% mutate( p_geno=lrt_pvalue(l_geno-l0,df=1),
                        p_interact=lrt_pvalue(l_interact-l_geno,df=df), 
                        p_joint=lrt_pvalue(l_interact-l0,df=df+1),
                        p_boot=lrt_pvalue(l_boot_interact - l_boot_geno, df ) ) %>%
    select(-starts_with("l"))
}

lrt_pvalue=function(halfdevi,df) pchisq(2*halfdevi,df=df,lower.tail = F)

bonferroni=function(g) { g %>% group_by(gene) %>% 
    summarize( cis_snp=cis_snp[which.min(p)], p=min(p) * length(p)  ) %>% 
    mutate( q=p %>% pmin(1) %>% p.adjust(method="BH") ) }

unscale=function(g) {
  sweep( sweep(g, 1, attr(g,"scaled:scale"), "*"), 1, attr(g,"scaled:center"), "+")
}


map_clusters_to_genes = function (intron_meta, exons_table) 
{
  foreach(chrom = sort(unique(intron_meta$chr)), .combine = bind_rows) %dopar% 
  {
    intron_chr = intron_meta %>% filter( chr == chrom )
    exons_chr = exons_table %>% filter( chr == chrom )
    three_prime_matches = inner_join(intron_chr, exons_chr, by=c(end="start") )
    five_prime_matches = inner_join(intron_chr, exons_chr, by=c(start="end") )
    all_matches = rbind( three_prime_matches %>% select( clu, gene_name ), 
                         five_prime_matches %>% select( clu, gene_name ) )  %>% 
      distinct()
  }
}

easy_impute=function(geno, prop_var=0.95) {
  temp=geno
  temp=t(scale(t(geno)))
  temp[is.na(temp)]=0
  s=svd(temp)
  v=s$d^2/sum(s$d^2)
  to_use=cumsum(v)<prop_var
  s$d[!to_use]=0.0
  recon=s$u %*% diag(s$d) %*% t(s$v)
  temp[is.na(geno)]=recon[is.na(geno)]
  temp=unscale(temp)
  stopifnot(max(abs(temp[!is.na(geno)]-geno[!is.na(geno)]))<1e-10)
  class(temp)="integer"
  temp
}

get_relatedness=function(filename, rna_inds) {
  ibd=read.table(filename, header=T)
  inds=intersect(ibd$Ind1, ibd$Ind2) # 1440!
  stopifnot(all(rna_inds %in% inds))
  ibd=ibd[ibd$Ind2 %in% rna_inds & ibd$Ind1 %in% rna_inds, ]
  ibd$Ind1=factor(ibd$Ind1, rna_inds)
  ibd$Ind2=factor(ibd$Ind2, rna_inds)
  
  errorCovariance=sparseMatrix(i=as.numeric(ibd$Ind1), j=as.numeric(ibd$Ind2), x=ibd$X0)
  #heatmap(as.matrix(errorCovariance), Rowv = NA, Colv=NA)
  #xor( t(errorCovariance)>0 , errorCovariance>0 ) # good
  errorCovariance=errorCovariance + t(errorCovariance) - diag(diag(errorCovariance))
  dimnames(errorCovariance)=list(rna_inds,rna_inds)
  det(errorCovariance) # 0.29 good
  as.matrix(errorCovariance)
}

require(irlba)
remove_PCs=function(input,num_PCs_to_remove) {
  # Remove PCs
  if (num_PCs_to_remove>0) {
    pca=irlba(as.matrix(input), nv=num_PCs_to_remove) # note this will remove dox signal
    recon=pca$u %*% (if (num_PCs_to_remove==1) pca$d else diag(pca$d)) %*% t(pca$v)
    #1 - mean((recon-input)^2) / mean(input^2) # 98% of variance is in first 5 PCs
    input=input-recon
  }
  input
}

require(doMC)

my_qqnorm=function(x) {
  n=length(x)
  a=if (n<=10) (3.0/8.0) else 0.5
  qnorm( (rank(x)-a)/(n+1.0-2.0*a) )
}

quantile_normalize_cols=function(input) {
  apply(input, 2, my_qqnorm)
}

quantile_normalize=function(input) {
  # Quantile normalization
  input_t=as.data.frame(t(input))
  res= foreach(l=as.list(input_t)) %dopar% {
    qqnorm(l,plot.it = F)$x
  }
  qnorm_input=t(as.matrix(as.data.frame(res)))
  dimnames(qnorm_input)=dimnames(input)
  qnorm_input
}

pvalue_qqplot=function(pvalues, nl10_obs_p_threshold=0) {
  n=length(pvalues)
  data.frame(x=-log10(seq(1/n, 1, length.out = n)), y=-log10(sort(pvalues)) ) %>%
    filter( y > nl10_obs_p_threshold ) %>%
    ggplot(aes(x,y)) + geom_point() + geom_abline(intercept=0,slope=1) + xlab("Expected -log10(p)") + ylab("Observed -log10(p)") + expand_limits(x=0, y=0)
}
