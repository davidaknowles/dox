// EAGLE extension for repeated measurements of the same individual, in T conditions. 
// BB(n,a,sigma(phi*beta*x)) is the model, where phi in (-1,1) corresponding to the phase. 
// This version allows for K SNPs
// No constraints on cofficients.
functions {
  // would be more efficient to pre-calc p
  real beta_binomial_reparam_lpmf(int y, int n, real g, real conc) {
    real p; 
    p = inv_logit(g);
    return beta_binomial_lpmf(y | n, conc*p, conc*(1.0-p));
  }
}
data {
  int<lower=0> N; // individuals
  int<lower=0> P; // covariates
  matrix[N,P] x; // covariates
  int<lower=0> ys[N]; // minor allele counts
  int<lower=0> ns[N]; // total coverage counts
  real<lower=0> concShape; // gamma prior on concentration parameters
  real<lower=0> concRate;
}
parameters {
  real<lower=0> conc; // concentration parameter
  vector[P] beta; // effect sizes
}
model {
  vector[N] xb; 
  xb = x * beta;
  for (n in 1:N) {
    target += beta_binomial_reparam_lpmf(ys[n] | ns[n], -xb[n], conc);
  }
  conc ~ gamma(concShape, concRate);
}
