data {
  int<lower=0> N; 
  vector[N] pvalues; 
  vector[N] weights;
}
parameters {
  real<lower=0,upper=1> pi0; // concentration parameter
  real<lower=0, upper=1> alpha; 
  real<lower=1> beta; 
}
model {
  alpha ~ gamma(1.001, 0.001);
  beta ~ gamma(1.001, 0.001);
  for(i in 1:N) {
    target += weights[i] * log_sum_exp(log(pi0) + beta_lpdf(pvalues[i] | 1, 1), log1m(pi0) + beta_lpdf(pvalues[i] | alpha, beta));
  }
}

