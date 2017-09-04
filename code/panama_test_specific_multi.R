

library("dplyr")
library("tidyr")
library(data.table)
source("utils.R")
#registerDoMC(7)
require(rstan)


source("load_data.R")

geneloc=read.table(paste0(DATADIR,"genelocGRCh38.txt.gz"),header=T,stringsAsFactors = F)

load("~/Dropbox/enviro_code/smalle_data/common_snps.RData")
hg38_snps = common_snps[,1:3]
rm(common_snps)
gc()
colnames(hg38_snps)=c("Ch","BP","RSID")

snploc=snploc %>% 
  mutate(Ch=substr(chr,4,nchar(chr)) %>% as.integer()) %>% 
  left_join(hg38_snps, by=c(Ch="Ch",pos="BP"))

#gwas_hit=data.frame(chr="chr15", pos=23457305, RSID="rs11855704")
normalization_approach="qq"
permuted="boot"

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

rownames(geneloc)=geneloc$geneid
cisdist=1e6

#errorhandling=if (interactive()) 'stop' else 'remove'
errorhandling='remove'

K=readRDS("../data/Kern.rds")
eig_K=eigen(K)

anno = anno %>% mutate(conc=as.factor(conc))
no_geno = model.matrix( ~ conc, data=anno) # [,2:5]

N=ncol(input)

mega_res=foreach(rsid=gwas_ld_filtered$other_snp, .combine = bind_rows) %dopar% {
  gwas_hit=snploc %>% filter(RSID==rsid)
  if (nrow(gwas_hit) != 1) return(NULL)

  cis_snp=snploc %>% filter(chr==gwas_hit$chr, pos==gwas_hit$pos) %>% .$snpid %>% as.character()

  genes = geneloc %>% filter(chr==gwas_hit$chr, gwas_hit$pos > (left-cisdist),  gwas_hit$pos < (right+cisdist) ) %>% .$geneid
  genes = intersect(rownames(input),genes)
  
  geno=genotype[cis_snp,anno$findiv]
  
  if (any(is.na(geno))) {
    snps=  snploc %>% filter(chr==gwas_hit$chr , abs(pos - gwas_hit$pos)<5e5) %>% .$snpid %>% as.character()
    imp_geno=easy_impute( genotype[snps,]  )
    geno=imp_geno[cis_snp,anno$findiv]
    #imp_geno[cis_snp,is.na(genotype[cis_snp,])]
  }
  
  foreach(gene=genes, .errorhandling=errorhandling, .combine = bind_rows) %dopar% {
    
    print(gene)
    y=input[gene,] %>% as.numeric
    y=y-mean(y)
    
    data=list(N=N,U_transpose_x=t(eig_K$vectors) %*% no_geno,P=ncol(no_geno), U_transpose_y=t(eig_K$vectors) %*% y %>% as.numeric, lambda=eig_K$values)
  
    init=list(sigma2=0.1, sigma2_k=1.0, beta=lm(y ~ no_geno - 1) %>% coef )
  
    fit_no_geno=optimizing(panama_test, data, init=init, as_vector=F)
    
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

df=4
mega_res=mega_res %>% mutate( p_geno=lrt_pvalue(l_geno-l0,df=1),
                      p_interact=lrt_pvalue(l_interact-l_geno,df=df), 
                      p_joint=lrt_pvalue(l_interact-l0,df=df+1),
                      p_boot=lrt_pvalue(l_boot_interact - l_boot_geno, df ) )

write.table(mega_res, file=paste0("../results/mega.txt"),sep="\t",row.names=F, quote=F)

hits = mega_res %>% group_by(cis_snp) %>% 
  top_n(1, -p_joint) %>% 
  ungroup() %>%
  mutate(cis_snp=as.integer(cis_snp)) %>% 
  left_join(snploc, by=c(cis_snp="snpid"))

mega_res %>% group_by(cis_snp) %>% summarize( p=min(p_joint) * length(p_joint), gene=gene[which.min(p_joint)] )
