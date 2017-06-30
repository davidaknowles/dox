
data {
  int<lower=0> N; # individuals
  int<lower=0> P; # covariates
  # U areeigenvectors of K 
  matrix[N,P] U_transpose_x; # rotated covariates
  vector[N] lambda; # eigenvalues of K
  vector[N] U_transpose_y;  # rotated expression
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
  R = U_transpose_y - U_transpose_x * beta; 
  target += - 0.5 * sum(R .* R ./ o );
  target += - 0.5 * sum(log(o));
}
