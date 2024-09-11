functions{
  real beta_neg_binomial_2_lpmf(array[] int y,vector mu, real sigma,real gamma);
  real beta_neg_binomial_2_lpmf(int y, real mu, real sigma, real gamma);
  int beta_neg_binomial_2_rng(real mu, real sigma, real gamma);
}
data {
  int<lower=0> N;           // number of data points
  int<lower=0> P;           // number of covariates
  matrix[N,P] X;            // covariates
  array[N] int<lower=0> y;  // target
  vector[N] offsett;        // offset (offset variable name is reserved)
}
parameters {
  real alpha;
  vector[P] beta;
  real<lower=0> sigma;
  real<lower=0> gamma;
}
model {
  // priors
  alpha ~ normal(0, 3);
  beta ~ normal(0, 3);
  gamma ~ normal(0, 3);
  sigma ~ normal(0, 3);
  // observation model
  vector[N] mu = offsett + exp(alpha + X * beta);
  target += beta_neg_binomial_2_lpmf(y| mu, sigma, gamma);
}
generated quantities {
  // log_lik for PSIS-LOO
  vector[N] log_lik;
  vector[N] y_rep;
  {
    vector[N] mu = offsett + exp(alpha + X * beta);
    for (n in 1:N) {

      log_lik[n] = beta_neg_binomial_2_lpmf(y[n]| mu[n], sigma, gamma);
      y_rep[n] = beta_neg_binomial_2_rng(mu[n], sigma, gamma);
    }
  }
}
