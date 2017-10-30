
data {
  int<lower=0> N; // individuals
  // U areeigenvectors of K 
  vector[N] lambda; // eigenvalues of K
  vector[N] U_transpose_y;  // rotated expression
}
parameters {
  real<lower=0> sigma2_k;
  real<lower=0> sigma2; 
}
model {
  vector[N] o; 
  vector[N] R; 
  o = sigma2_k * lambda + sigma2; 
  target += - 0.5 * sum(U_transpose_y .* U_transpose_y ./ o );
  target += - 0.5 * sum(log(o));
}
