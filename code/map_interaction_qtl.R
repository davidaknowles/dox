require(rstan)

## Wrapper to supress extraneous stdout from rstan::optimizing
stan_optimizing_wrapper = function(...) {
  zz <- file( "/dev/null", open = "wt")
  sink(zz)
  sink(zz, type = "message")
  res=optimizing(...)
  sink(type="message")
  sink()
  res
}

map_interaction_qtl = function(input, genotype, geneloc, snploc, anno, sample_kernel, normalization_approach, permutation_approach, cisdist, checkpoint_dir) {
  
  #input=remove_PCs(input, num_PCs_to_remove)
  if (normalization_approach=="qq") {
    input=quantile_normalize(input)
  } else if (normalization_approach=="log") {
    input=input %>% t %>% scale %>% t
  } else if (normalization_approach=="linear") {
    input=(2^input) %>% t %>% scale %>% t
  } else if (normalization_approach !="none" ) {
    stop(paste0("Invalid normalization_approach:",normalization_approach))
  }
  
  genes=intersect(rownames(input),geneloc$geneid)
  rownames(geneloc)=geneloc$geneid
  
  eigen_sample_kernel=eigen(sample_kernel)
  
  anno = anno %>% mutate(conc=as.factor(conc))
  no_geno = model.matrix( ~ conc, data=anno) # [,2:5]
  
  N=ncol(input)
  
  panama_test=stan_model("panama_test.stan")

  errorhandling=if (interactive()) 'stop' else 'remove'
  
  foreach(gene=genes, .errorhandling=errorhandling, .combine = bind_rows) %do% {
    
    check_fn=paste0(checkpoint_dir,gene,".txt.gz")
    if (file.exists(check_fn)) {
      return(read.table(check_fn, header=T, stringsAsFactors = F, sep="\t"))
    }
    
    cis_snps=snploc %>% 
      filter(chr==geneloc[gene,"chr"], (geneloc[gene,"left"]-cisdist) < pos,  (geneloc[gene,"right"]+cisdist) > pos) %>%
      .$snpid %>%
      as.character()
    cat(gene,length(cis_snps)," cis snps\n")
    
    if (length(cis_snps)==0) return(NULL)
    print(1)
    y=input[gene,] %>% as.numeric
    y=y-mean(y)
    
    imp_geno=easy_impute(genotype[cis_snps,,drop=F])
    print(2)
    if (permutation_approach=="permute") colnames(imp_geno)=colnames(imp_geno)[ sample(ncol(imp_geno),ncol(imp_geno)) ]
    # cis_snp=as.character(cis_snps)[1]
    print("2a")    
    data=list(N=N,U_transpose_x=t(eigen_sample_kernel$vectors) %*% no_geno,P=ncol(no_geno), U_transpose_y=t(eigen_sample_kernel$vectors) %*% y %>% as.numeric, lambda=eigen_sample_kernel$values)
    print("2b")
    init=list(sigma2=0.1, sigma2_k=1.0, beta=lm(y ~ no_geno - 1) %>% coef )
    print(3)
    fit_no_geno=stan_optimizing_wrapper(panama_test, data, init=init, as_vector=F)
    print("Here")
    gene_results = foreach(cis_snp=cis_snps, .errorhandling=errorhandling, .combine = bind_rows) %do% {
    	print(cis_snp)
      geno=imp_geno[cis_snp,anno$findiv]
      if (sum(imp_geno[cis_snp,]) < 5.0) {
      	 print("Too few minor alleles")
	 return(NULL)
 	 }
      
      lrt = function(data) {
        data$U_transpose_x=t(eigen_sample_kernel$vectors) %*% cbind( no_geno, geno )
        data$P=ncol(data$U_transpose_x)
        init=fit_no_geno$par
        init$beta=c(init$beta,0.0)
        
        fit_geno=stan_optimizing_wrapper(panama_test, data, init=init, as_vector=F )
        
        interact=model.matrix(~geno:conc,data=anno)
        interact=interact[,3:ncol(interact)]
        data_interact=data
        data_interact$U_transpose_x=t(eigen_sample_kernel$vectors) %*% cbind( no_geno, geno, interact )
        data_interact$P=ncol(data_interact$U_transpose_x)
        
        init=fit_geno$par
        init$beta=c(init$beta,numeric(ncol(interact)))
        fit_interact=stan_optimizing_wrapper(panama_test, data_interact, init=init, as_vector=F)
        
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
  print("Saving results")
  checkpoint_file= gzfile( check_fn,"w")
  gene_results %>% format(digits=5) %>% write.table(checkpoint_file, quote = F, row.names = F, col.names = T, sep="\t")
  close(checkpoint_file)
  
  gene_results
  }
}

# Could just zcat everything into one file rather than bind_rows

