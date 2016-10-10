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
model {
  vector[N] means[K];
  for (k in 1:K)
    means[k] <- x * betas[k];
  p ~ dirichlet( rep_vector(1.0/K,K) ); 
  for (d in 1:D) {
    vector[K] logprob; 
    for (k in 1:K) {
      real l[N]; 
      for (n in 1:N)
        l[n] <- normal_log(y[d,n], means[k][n], sigma);
      logprob[k] <- sum(l);
    }
    increment_log_prob(log_sum_exp(log(p)+logprob)); 
  }
}
