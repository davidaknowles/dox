
data {
  int<lower=0> N; # individuals
  int<lower=0> P; # covariates
  int<lower=0> Pfix; # covariates
  matrix[N,N] x[P]; 
  matrix[N,Pfix] xfix; 
  vector[N] y; 
}
parameters {
  real<lower=0> s[P]; 
  vector[Pfix] beta; 
}
model {
  matrix[N,N] sigma;
  sigma <- rep_matrix(0.0, N, N);
  for (p in 1:P) {
    sigma <- sigma + s[p] * x[p]; 
  }
  y ~ multi_normal( xfix * beta, sigma);
}
