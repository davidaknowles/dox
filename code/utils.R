require(Matrix)

unscale=function(g) {
  sweep( sweep(g, 1, attr(g,"scaled:scale"), "*"), 1, attr(g,"scaled:center"), "+")
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
quantile_normalize=function(input) {
  # Quantile normalization
  input_t=as.data.frame(t(input))
  res= foreach(l=as.list(input_t)) %do% {
    qqnorm(l,plot.it = F)$x
  }
  qnorm_input=t(as.matrix(as.data.frame(res)))
  dimnames(qnorm_input)=dimnames(input)
  qnorm_input
}