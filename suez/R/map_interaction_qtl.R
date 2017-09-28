
## Wrapper to supress extraneous stdout from rstan::optimizing
## Doing this seems to open too many file handles though :(
stan_optimizing_wrapper = function(...) {
  #zz <- file( "/dev/null", open = "wt")
  #sink(zz)
  #sink(zz, type = "message")
  res=optimizing(...)
  #sink(type="message")
  #sink()
  res
}

#' Simple SVD based imputation of missing genotypes
#' 
#' @param geno [samples] x [SNPs] genotype matrix (0/1/2)
#' @param prop_var Proportion of variance that the PCs should explain
#' 
#' @return Complete genotype matrix. 
easy_impute=function(geno, prop_var=0.95) {
  temp=geno
  temp=t(scale(t(geno)))
  temp[is.na(temp)]=0
  s=svd(temp)
  v=s$d^2/sum(s$d^2)
  to_use=cumsum(v)<prop_var
  s$d[!to_use]=0.0
  recon=s$u %*% fix_diag(s$d) %*% t(s$v)
  temp[is.na(geno)]=recon[is.na(geno)]
  temp=unscale(temp)
  stopifnot(max(abs(temp[!is.na(geno)]-geno[!is.na(geno)]))<1e-10)
  temp=round(temp)
  class(temp)="integer"
  temp
}

#' qqnorm without plotting
#' 
#' @param x Numeric vector input
#' @return Quantile normalized version
qqnorm_no_plot=function(x) {
  n=length(x)
  a=if (n<=10) (3.0/8.0) else 0.5
  qnorm( (rank(x)-a)/(n+1.0-2.0*a) )
}

#' Quantile normalize columns
quantile_normalize_cols=function(input) {
  apply(input, 2, qqnorm_no_plot)
}

#' Quantile normalize rows (to normal)
#' @import magrittr
quantile_normalize=function(input) {
  apply(input, 1, qqnorm_no_plot) %>% t()
}

#' Response (e)QTL mapping
#' 
#' @param input [genes x samples] matrix of log_2 expression values, e.g. cpm or fpkm
#' @param genotype [SNPs x individuals] matrix of genotypes (can have some NAs)
#' @param geneloc [genes x 4] data.frame with columns (geneid, chr, left, right). The intersect of geneloc$geneid and rownames(input) will be tested. 
#' @param snploc [SNPs x 3] data.frame with columns (snpid, chr, pos) where snpid's map to rows of `genotype`
#' @param anno [samples x 2] data.frame with columns (individual, condition) where individual must correspond to column names of `genotype`
#' @param sample_kernel [samples x samples] covariance matrix learned in suez step 1. If you don't want to control for latent confounders leave this as NULL to use the matrix representing which samples are from the same individual. 
#' @param normalization_approach One of (qq,log,linear) determining how to normalize within each gene. We recommend qq (quantile normalization to normal) to ensure the assumptions of the test hold.
#' @param permutation_approach One of ("none","permute","boot"). We recommend boot (parametric bootstrap) since the permutation test is not really valid for interaction effect testing. 
#' @param cisdist How far from the gene boundaries to look for SNPs to test
#' @param checkpoint_dir An optional directory to store per gene results. If mapping crashes these results will be reused to save time. 
#' @param debug Controls whether you see stan::optimizing messages and whether errors are suppressed. 
#' @importFrom rstan optimizing
#' @import dplyr
#' @import doMC
#' @useDynLib suez, .registration = TRUE 
#' @export
map_interaction_qtl = function(input, genotype, geneloc, snploc, anno, sample_kernel=NULL, normalization_approach="qq", permutation_approach="boot", cisdist=1e5, checkpoint_dir=NULL, debug=F) {
  
  if (normalization_approach=="qq") {
    input=quantile_normalize(input)
  } else if (normalization_approach=="log") {
    input=input %>% t %>% scale %>% t
  } else if (normalization_approach=="linear") {
    input=(2^input) %>% t %>% scale %>% t
  } else if (normalization_approach !="none" ) {
    stop(paste0("Invalid normalization_approach:",normalization_approach))
  }
  
  if (!is.null(checkpoint_dir)) dir.create(checkpoint_dir, recursive = T, showWarnings = F)
  
  genes=intersect(rownames(input),geneloc$geneid)
  rownames(geneloc)=geneloc$geneid
  
  if (is.null(sample_kernel)) {
    sample_kernel=outer(anno$individual,anno$individual,"==")
    class(sample_kernel)="numeric"
  }
  
  eigen_sample_kernel=eigen(sample_kernel)
  
  anno = anno %>% mutate(condition=as.factor(condition))
  no_geno = model.matrix( ~ condition, data=anno) # [,2:5]
  
  N=ncol(input)
  
  errorhandling=if (debug) 'stop' else 'remove'
  
  if (!debug) { 
    zz <- file( "/dev/null", open = "wt")
    sink(zz)
    sink(zz, type = "message")
    on.exit({
      sink(type="message")
      sink()
    })
  }
  
  foreach(gene=genes, .errorhandling=errorhandling, .combine = bind_rows) %do% {
    
    if (!is.null(checkpoint_dir)){
      check_fn=paste0(checkpoint_dir,gene,".txt.gz")
      if (file.exists(check_fn)) {
        return(read.table(check_fn, header=T, stringsAsFactors = F, sep="\t"))
      }
    }
    
    cis_snps=snploc %>% 
      filter(chr==geneloc[gene,"chr"], 
             (geneloc[gene,"left"]-cisdist) < pos,  
             (geneloc[gene,"right"]+cisdist) > pos) %>%
      .$snpid %>%
      as.character()
    cat(gene,length(cis_snps)," cis snps\n")
    
    if (length(cis_snps)==0) return(NULL)

    y=input[gene,] %>% as.numeric
    y=y-mean(y)
    
    imp_geno=easy_impute(genotype[cis_snps,,drop=F])

    if (permutation_approach=="permute") colnames(imp_geno)=colnames(imp_geno)[ sample(ncol(imp_geno),ncol(imp_geno)) ]

    data=list(N=N,U_transpose_x=t(eigen_sample_kernel$vectors) %*% no_geno,P=ncol(no_geno), U_transpose_y=t(eigen_sample_kernel$vectors) %*% y %>% as.numeric, lambda=eigen_sample_kernel$values)

    init=list(sigma2=0.1, sigma2_k=1.0, beta=lm(y ~ no_geno - 1) %>% coef )

    fit_no_geno=stan_optimizing_wrapper(stanmodels$suez_step_2, data, init=init, as_vector=F)

    gene_results = foreach(cis_snp=cis_snps, .errorhandling=errorhandling, .combine = bind_rows) %dopar% {

      geno=imp_geno[cis_snp,anno$individual]
      if (sum(imp_geno[cis_snp,]) < 5.0) { return(NULL) }

      lrt = function(data) {
        data$U_transpose_x=t(eigen_sample_kernel$vectors) %*% cbind( no_geno, geno )
        data$P=ncol(data$U_transpose_x)
        init=fit_no_geno$par
        init$beta=c(init$beta,0.0)
        
        fit_geno=stan_optimizing_wrapper(stanmodels$suez_step_2, data, init=init, as_vector=F )
        
        interact=model.matrix(~geno:condition,data=anno)
        interact=interact[,3:ncol(interact)]
        data_interact=data
        data_interact$U_transpose_x=t(eigen_sample_kernel$vectors) %*% cbind( no_geno, geno, interact )
        data_interact$P=ncol(data_interact$U_transpose_x)
        
        init=fit_geno$par
        init$beta=c(init$beta,numeric(ncol(interact)))
        fit_interact=stan_optimizing_wrapper(stanmodels$suez_step_2, data_interact, init=init, as_vector=F)
        
        list( fit_geno=fit_geno, fit_interact=fit_interact )
      }
      
      lrt_true=lrt( data )
      fit_geno=lrt_true$fit_geno
      fit_interact=lrt_true$fit_interact
      
      res=data.frame(gene=gene, cis_snp=cis_snp, l0=fit_no_geno$value, l_geno=fit_geno$value, l_interact=fit_interact$value, stringsAsFactors=F  )
      
      if (permutation_approach=="boot") {
        Sigma = fit_geno$par$sigma2_k * sample_kernel + fit_geno$par$sigma2 * diag(N)
        chol_Sigma = chol(Sigma)
        xb = cbind( no_geno, geno ) %*% fit_geno$par$beta
        y_boot = t(chol_Sigma) %*% rnorm(N) + xb
        data_boot = data
        data_boot$U_transpose_y = t(eigen_sample_kernel$vectors) %*% y_boot %>% as.numeric()
        lrt_boot = lrt( data_boot )
        res$l_boot_geno=lrt_boot$fit_geno$value
        res$l_boot_interact=lrt_boot$fit_interact$value
      }
      
      res
    }
    if (!is.null(checkpoint_dir)){
      print("Saving results")
      checkpoint_file= gzfile( check_fn,"w")
      gene_results %>% format(digits=5) %>% write.table(checkpoint_file, quote = F, row.names = F, col.names = T, sep="\t")
      close(checkpoint_file)
    }
    
    gene_results
  }
}

# Could just zcat everything into one file rather than bind_rows

