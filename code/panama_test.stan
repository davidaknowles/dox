
data {
  int<lower=0> N; # individuals
  int<lower=0> P; # covariates
  matrix[N,N] U_transpose; # eigenvectors of K 
  vector[N] lambda; # eigenvalues of K
  vector[N] y; 
  matrix[N,P] x;
}
parameters {
  real<lower=0> sigma2_k;
  real<lower=0> sigma2; 
  vector[P] beta; 
}
model {
  vector[N] o; 
  vector[N] R; 
  o = sigma2_k * lambda + sigma2; 
  R = U_transpose * (y - x * beta); 
  target += - 0.5 * sum(R .* R ./ o );
  target += - 0.5 * sum(log(o));
}
