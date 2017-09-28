
data {
  int<lower=0> N; // individuals
  int<lower=0> G; // genes
  int<lower=0> P; // covariates
  int<lower=0> P_adjustable; // covariates
  matrix[N,N] x[P]; 
  vector[N] y[G]; 
}
parameters {
  real<lower=0> s[P];
  real<lower=0> s_adjustable[P_adjustable]; 
  unit_vector[N] x_adjustable[P_adjustable]; 
}
model {
  matrix[N,N] sigma;
  matrix[N,N] L; 
  sigma = rep_matrix(0.0, N, N);
  for (p in 1:P) 
    sigma = sigma + s[p] * x[p]; 
  for (p in 1:P_adjustable) 
    sigma = sigma + s_adjustable[p] * ( x_adjustable[p] * to_row_vector(x_adjustable[p]) ); 
  L=cholesky_decompose(sigma);
  for (g in 1:G)
    y[g] ~ multi_normal_cholesky( rep_vector(0,N), L);
}
