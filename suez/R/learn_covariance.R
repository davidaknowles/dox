

#' Function to learn covariance structure across samples using latent factors
#' 
#' @param input [genes x samples] Normalized (e.g. quantile normalized) matrix of log2 gene expression
#' @param anno [samples x 2] data.frame with `individual` and `condition` columns. If a kinship matrix is provided the row and column names must correspond to the values used in the `individual` column. 
#' @param P Number of latent factors to model. suez can prune out unnecessary factors (in principle). 
#' @param kinship_matrix Representing relatedness between individuals. Optional. 
#' @param iterations Max iterations of LBFGS to run
#' @param same_ind_initial_weight Initialization. 
#' @param kinship_matrix_initial_weight Initialization.
#' @import rstan
#' @export
learn_covariance = function(input, anno, kinship_matrix=NULL, P=20, iterations=200, same_ind_initial_weight=0.01, kinship_matrix_initial_weight=0.01 )
{
  
  self_outer=function(g) outer(g,g)
  
  # Regress out condition effect
  X = model.matrix( ~ as.factor(condition), data=anno)
  b=solve( t(X) %*% X, t(X) %*% t(input) )
  resi = input - t(X %*% b)
  
  svd_ge=svd(input)
  
  # Latent factors initialized using SVD
  x_adjustable=foreach(p=1:P) %do% { svd_ge$v[,p] }
  sinit=svd_ge$d^2 / nrow(input)
  
  same_ind=outer(anno$individual, anno$individual, "==") * 1
  same_condition=outer(anno$condition, anno$condition, "==") * 1
  
  variance_components=list( isotropic_noise=diag(ncol(input)), same_ind=same_ind, same_condition=same_condition)
  # residual after removing first P PCs (used to estimate noise variance)
  rr = resi - svd_ge$u[,1:P] %*% diag(svd_ge$d[1:P]) %*% t(svd_ge$v[,1:P])
  variance_component_weights=c( mean(rr^2), same_ind_initial_weight, mean(b^2) )
  
  if (!is.null(kinship_matrix)) {
    variance_components$kinship_matrix=kinship_matrix[ anno$individual, anno$individual ]
    variance_component_weights=c(variance_component_weights, kinship_matrix_initial_weight)
  }

  dat=list(y=resi, N=ncol(input), G=nrow(input), P=length(variance_components), x=variance_components, P_adjustable=length(x_adjustable))

  init=list( s_adjustable=sinit[1:P], s=variance_component_weights, x_adjustable=x_adjustable )
  
  optimized_fit=optimizing(stanmodels$suez_step_1, data=dat, init=init, verbose=T, as_vector=F, iter=iterations)
  
  names(optimized_fit$par$s) = names(variance_components)
  
  cov_suez = Reduce("+", foreach(p=1:dat$P_adjustable) %do% { 
    optimized_fit$par$s_adjustable[p] * self_outer(optimized_fit$par$x_adjustable[p,]) 
    } )
  
  cov_suez = cov_suez + optimized_fit$par$s["same_ind"] * variance_components[["same_ind"]]
  if (!is.null(kinship_matrix)) {
    cov_suez = cov_suez + optimized_fit$par$s["kinship_matrix"] * variance_components[["kinship_matrix"]]
  }
  
  cov_suez
}