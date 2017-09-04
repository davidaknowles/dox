
library("dplyr")
library("tidyr")
library(data.table)
source("utils.R")
registerDoMC(16)
source("load_data.R")
require(rstan)

if (interactive()) {
  chrom="chr15"
  normalization_approach="qq"
  permuted="boot"
  cisdist=1e5
} else {
  ca=commandArgs(trailingOnly = T)
  chrom=ca[2]
  permuted=ca[3]
  normalization_approach=ca[1]
  cisdist=as.numeric(ca[4])
}

geneloc=geneloc[geneloc$chr==chrom,]
snploc=snploc[snploc$chr==chrom,]

stopifnot(all(as.character(snploc$snpid) %in% rownames(genotype) ))
genotype=genotype[as.character(snploc$snpid),]

#input=remove_PCs(input, num_PCs_to_remove)
if (normalization_approach=="qq") {
  input=quantile_normalize(input)
} else if (normalization_approach=="log") {
  input=input %>% t %>% scale %>% t
} else if (normalization_approach=="linear") {
  input=(2^input) %>% t %>% scale %>% t
}

panama_test=stan_model("panama_test.stan")

genes=intersect(rownames(input),geneloc$geneid)
rownames(geneloc)=geneloc$geneid
errorhandling=if (interactive()) 'stop' else 'remove'

K=readRDS("../data/Kern.rds")
eig_K=eigen(K)

anno = anno %>% mutate(conc=as.factor(conc))
no_geno = model.matrix( ~ conc, data=anno) # [,2:5]

N=ncol(input)

results=foreach(gene=genes, .errorhandling=errorhandling, .combine = bind_rows) %do% {
  
  print(gene)
  y=input[gene,] %>% as.numeric
  y=y-mean(y)
  
  cis_snps=snploc[ ((geneloc[gene,"left"]-cisdist) < snploc$pos) & ((geneloc[gene,"right"]+cisdist) > snploc$pos), "snpid" ]
  cis_snps=as.character(cis_snps)

  imp_geno=easy_impute(genotype[cis_snps,])
  
  if (permuted=="permute") colnames(imp_geno)=colnames(imp_geno)[ sample(ncol(imp_geno),ncol(imp_geno)) ]
  # cis_snp=as.character(cis_snps)[1]
  
  data=list(N=N,U_transpose_x=t(eig_K$vectors) %*% no_geno,P=ncol(no_geno), U_transpose_y=t(eig_K$vectors) %*% y %>% as.numeric, lambda=eig_K$values)

  init=list(sigma2=0.1, sigma2_k=1.0, beta=lm(y ~ no_geno - 1) %>% coef )

  fit_no_geno=optimizing(panama_test, data, init=init, as_vector=F)
  foreach(cis_snp=cis_snps, .errorhandling=errorhandling, .combine = bind_rows) %dopar% {
    geno=imp_geno[cis_snp,anno$findiv]
    if (sum(imp_geno[cis_snp,]) < 5.0) return(NULL)

    lrt = function(data) {
      data$U_transpose_x=t(eig_K$vectors) %*% cbind( no_geno, geno )
      data$P=ncol(data$U_transpose_x)
      init=fit_no_geno$par
      init$beta=c(init$beta,0.0)
      
      fit_geno=optimizing(panama_test, data, init=init, as_vector=F )
      
      interact=model.matrix(~geno:conc,data=anno)
      interact=interact[,3:ncol(interact)]
      data_interact=data
      data_interact$U_transpose_x=t(eig_K$vectors) %*% cbind( no_geno, geno, interact )
      data_interact$P=ncol(data_interact$U_transpose_x)
    
      init=fit_geno$par
      init$beta=c(init$beta,numeric(ncol(interact)))
      fit_interact=optimizing(panama_test, data_interact, init=init, as_vector=F)
      
      list( fit_geno=fit_geno, fit_interact=fit_interact )
    }
    
    lrt_true=lrt( data )
    fit_geno=lrt_true$fit_geno
    fit_interact=lrt_true$fit_interact
    
    res=data.frame(gene=gene, cis_snp=cis_snp, l0=fit_no_geno$value, l_geno=fit_geno$value, l_interact=fit_interact$value  )
    
    if (permuted=="boot") {
      Sigma = fit_geno$par$sigma2_k * K + fit_geno$par$sigma2 * diag(N)
      chol_Sigma = chol(Sigma)
      xb = cbind( no_geno, geno ) %*% fit_geno$par$beta
      y_boot = t(chol_Sigma) %*% rnorm(N) + xb
      data_boot = data
      data_boot$U_transpose_y = t(eig_K$vectors) %*% y_boot %>% as.numeric()
      lrt_boot = lrt( data_boot )
      res$l_boot_geno=lrt_boot$fit_geno$value
      res$l_boot_interact=lrt_boot$fit_interact$value
    }
    
    #lrt_interact=2.0*(fit_interact$value - fit_geno$value)
    #df_interact=ncol(interact)
    # pchisq(lrt,df,lower.tail=F)
    
    res
  }
  
}

resdir=paste0(DATADIR,"panama_",normalization_approach,"_",permuted,"_",cisdist,"/")
dir.create(resdir)

gz1 = gzfile(paste0(resdir,chrom,".txt.gz"),"w")
results %>% format(digits=5) %>% write.table(gz1, quote = F, row.names = F, col.names = T, sep="\t")
close(gz1)
