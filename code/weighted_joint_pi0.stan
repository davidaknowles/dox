data {
  int<lower=0> N; 
  vector[N] pvalues;
  matrix[N,N] weights;
}
parameters {
  simplex[4] pi0; 
  real<lower=0, upper=1> alpha[2]; 
  real<lower=1> beta[2]; 
}
model {
  alpha ~ gamma(1.001, 0.001);
  beta ~ gamma(1.001, 0.001);
  for(i in 1:N) {
    for(j in 1:N) {
      vector[4] temp;
      //temp = log(pi0);
      temp[1] = log(pi0[1]) + beta_lpdf(pvalues[i] | 1, 1) +  beta_lpdf(pvalues[j] | 1, 1); 
      temp[2] = log(pi0[2]) + beta_lpdf(pvalues[i] | alpha[1], beta[1]) +  beta_lpdf(pvalues[j] | 1,1); 
      temp[3] = log(pi0[3]) + beta_lpdf(pvalues[i] | 1, 1) +  beta_lpdf(pvalues[j] | alpha[2], beta[2]); 
      temp[4] = log(pi0[4]) + beta_lpdf(pvalues[i] | alpha[1], beta[1]) +  beta_lpdf(pvalues[j] | alpha[2], beta[2]); 
      target += weights[i,j] * log_sum_exp(temp);
    }
  }
}
