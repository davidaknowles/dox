data {
  int<lower=0> N; # samples (217)
  int<lower=0> P; # covariates
  int<lower=0> D; # genes (1e4?)
  int<lower=0> K; 
  matrix[N,P] x; 
  matrix[D,N] y;
}
parameters {
  simplex[K] p; 
  real<lower=0> sigma; 
  vector[P] betas[K];
}
transformed parameters {
  vector[N] means[K];
  vector[K] logprob[D];
  vector[K] logp; 
  for (k in 1:K) {
    means[k] <- x * betas[k];
    logp[k] <- log(p[k]);
  }
  for (d in 1:D) {
    for (k in 1:K) {
      real l[N]; 
      for (n in 1:N)
        l[n] <- normal_log(y[d,n], means[k][n], sigma);
      logprob[d][k] <- sum(l) + logp[k];
    }
  }
}
model {
  p ~ dirichlet( rep_vector(1.0/K,K) ); 
  for (d in 1:D) {
    increment_log_prob(log_sum_exp(logprob[d])); 
  }
}
