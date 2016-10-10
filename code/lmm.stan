
data {
  int<lower=0> N; # individuals
  int<lower=0> P; # covariates
  matrix[N,N] x[P]; 
  vector[N] y; 
}
parameters {
  real<lower=0> s[P]; 
}
model {
  matrix[N,N] sigma;
  sigma <- rep_matrix(0.0, N, N);
  for (p in 1:P) {
    sigma <- sigma + s[p] * x[p]; 
  }
  y ~ multi_normal( rep_vector(0,N), sigma);
}
