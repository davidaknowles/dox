require(Matrix)

cbPalette <- c( "#009E73","#F0E442","#D55E00", "#999999", "#E69F00", "#56B4E9",  "#0072B2",  "#CC79A7")

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

fix_diag=function(x) {
  if(length(x)==1) matrix(x) else diag(x)
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
  apply(input, 1, my_qqnorm) %>% t()
}

pvalue_qqplot=function(pvalues, nl10_obs_p_threshold=0) {
  n=length(pvalues)
  data.frame(x=-log10(seq(1/n, 1, length.out = n)), y=-log10(sort(pvalues)) ) %>%
    filter( y > nl10_obs_p_threshold ) %>%
    ggplot(aes(x,y)) + geom_point() + geom_abline(intercept=0,slope=1) + xlab("Expected -log10(p)") + ylab("Observed -log10(p)") + expand_limits(x=0, y=0)
}

get_qqplot_data=function(pvalues) {
  n=length(pvalues)
  data.frame(x=-log10(seq(1/n, 1, length.out = n)), y=-log10(sort(pvalues)) )
}

pvalue_qqplot=function(pvalues, nl10_obs_p_threshold=0) {
  get_qqplot_data(pvalues) %>%
    filter( y > nl10_obs_p_threshold ) %>%
    ggplot(aes(x,y)) + geom_point() + geom_abline(intercept=0,slope=1) + xlab("Expected -log10(p)") + ylab("Observed -log10(p)") + expand_limits(x=0, y=0)
}

pvalue_qqplot_multi=function(pvalues, nl10_obs_p_threshold=0) {
  pvalues %>% group_by(group) %>% mutate(x=-log10(seq(1/n(), 1, length.out = n())), y=-log10(sort(p)) ) %>% ungroup() %>%
    filter( y > nl10_obs_p_threshold ) %>%
    ggplot(aes(x,y,col=group)) + geom_point() + geom_abline(intercept=0,slope=1) + xlab("Expected -log10(p)") + ylab("Observed -log10(p)") 
}

pvalue_qqplot_multi_thin=function(pvalues, nl10_obs_p_threshold=0, fidelity=1000) {
  pvalues %>% 
    group_by(group) %>% 
    mutate(x=-log10(seq(1/n(), 1, length.out = n())), y=-log10(sort(p)) ) %>% 
    ungroup() %>%
    filter( y > nl10_obs_p_threshold ) %>%
    mutate(x_r=round(x*fidelity), y_r=round(y*fidelity)) %>%
    group_by(group, x_r, y_r) %>% 
    slice(which.min(p))  %>%
    ungroup() %>% 
    ggplot(aes(x,y,col=group)) + geom_point() + geom_abline(intercept=0,slope=1) + xlab("Expected -log10(p)") + ylab("Observed -log10(p)") 
}


