my_pi0est = function (p, 
          lambda = seq(0.05, 0.95, 0.05), 
          smooth.df = 3) 
{
  if (length(p)==0) return(list(pi0=NA, pi0.lambda = NA*lambda, lambda = lambda))
  p <- p[!is.na(p)]
  lambda <- sort(lambda)
  ll <- length(lambda)
  pi0.lambda <- sapply(lambda, function(l) mean(p >= l)/(1 - l))
  spi0 <- smooth.spline(lambda, pi0.lambda, df = smooth.df)
  pi0Smooth <- predict(spi0, x = lambda)$y
  pi0 <- pi0Smooth[length(pi0Smooth)]
  return(list(pi0 = pi0, pi0.lambda = pi0.lambda, lambda = lambda))
}

my_joint_pi0est = function (p1, p2, 
                      lambda = seq(0.05, 0.95, 0.05), 
                      smooth.df = 3) 
{
  pi0_1 = my_pi0est(p1, lambda, smooth.df)$pi0
  pi0_2 = my_pi0est(p2, lambda, smooth.df)$pi0
  to_keep=(!is.na(p1)) & (!is.na(p2))
  p1 <- p1[to_keep]
  p2 <- p2[to_keep]
  lambda <- sort(lambda)
  pi0.lambda <- sapply(lambda, function(l) (mean( (p1 >= l) & (p2 >= l) ) /((1 - l)^2)))
  spi0 <- smooth.spline(lambda, pi0.lambda, df = smooth.df)
  pi0Smooth <- predict(spi0, x = lambda)$y
  pi0 <- pi0Smooth[length(pi0Smooth)]
  #  jaccard=(1-pi0) / ((1-pi0_1) + (1-pi0_2) - (1-pi0)))
  jaccard=( (1-pi0_1) + (1-pi0_2) - (1-pi0))  / (1-pi0)
  return(list(pi0 = pi0, pi0_1=pi0_1, pi0_2=pi0_2, pi0.lambda = pi0.lambda, lambda = lambda, jaccard=jaccard))
}

#joint_pi0_stan=stan_model("joint_pi0.stan")

#mix_pi0_stan=stan_model("mix_pi0.stan")

require(rstan)
require(qvalue)
require(gplots)

mix_weighted_pi0_stan=stan_model("mix_pi0.stan")
weighted_joint_pi0_stan=stan_model("weighted_joint_pi0.stan")

#temp=ase_qtls %>% head(30e6) %>% filter(!is.na(p_geno), !is.na(p_interact)) %>% mutate(p_geno=pmin(p_geno,.999), p_interact=pmin(p_interact,.999))


joint_pi0_estimator=function(pvalues, breaks=seq(0,1,by=0.01), restarts=5) {
  
  pvalues[,1:2]=pmin(pvalues[,1:2])
  
  init=foreach(pv=as.list(pvalues)[1:2]) %do% {
    pi0_init = pi0est(pv)$pi0
    h=hist(pv, breaks=breaks, plot=F)
    fits=foreach(i=seq_len(restarts)) %do% {
      h %$% optimizing(mix_weighted_pi0_stan, data=list(N=length(density), weights=density, pvalues=mids), init=list(pi0=pi0_init), as_vector=F)
    }
    fits[[ which.max(foreach(fit=fits) %do% { fit$value } ) ]]$par
  }
  
  pi0_init=c(init[[1]]$pi0*init[[2]]$pi0, (1-init[[1]]$pi0)*init[[2]]$pi0, init[[1]]$pi0*(1-init[[2]]$pi0), (1-init[[1]]$pi0)*(1-init[[2]]$pi0))
  
  N=length(breaks)-1
  h=gplots::hist2d(pvalues[,1], pvalues[,2], N)
  alpha=foreach(ini=init, .combine=c) %do% { ini$alpha }
  beta=foreach(ini=init, .combine=c) %do% { ini$beta }
  
  fits=foreach(i=seq_len(restarts)) %do% {
    h %$% optimizing(weighted_joint_pi0_stan, 
                     data=list(N=N, 
                               pvalues=x, 
                               weights=counts/sum(counts)), 
                     init=list(pi0=pi0_init, 
                               alpha=alpha, 
                               beta=beta), 
                     as_vector=F)
  }
  o=fits[[ which.max(foreach(fit=fits) %do% { fit$value } ) ]]
  
  #o$par$pi0[4] - (o$par$pi0[2]+o$par$pi0[4])*(o$par$pi0[3]+o$par$pi0[4]) # covariance
  list(jaccard=o$par$pi0[4] / sum(o$par$pi0[2:4]),  # jaccard like index
       fit=o) 
}

