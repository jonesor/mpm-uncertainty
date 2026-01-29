
data {
  int<lower=0> N; // number of data points, length(y)
  int<lower=0> K; // number of monthly lags
  int y[N];
  vector[N] offset;
  matrix[N,K] X;  // matrix of climate covariates
}

transformed data {
  real month[K]; // vector of months, to create distance matrix
  for(k in 1:K) month[k] = k;
}

parameters {
  real<lower=0> eta; // maximum covariance for betas
  real<lower=0> rho; // degree of temporal autocorrelation for betas
  real mu_alpha;
  real<lower=0> sigma_alpha;
  real mu_beta;
  vector[N] z1;      // unit normal prior for non-centered alpha term
  vector[K] z2;      // unit normal prior for non-centered beta term
}

transformed parameters {
  vector[N] lambda;
  vector[N] alpha_year;
  vector[K] beta;
  matrix[K,K] Sigma; // covariance matrix
  matrix[K,K] L;     // cholesky of covariance matrix
  
  // covariance
  Sigma = cov_exp_quad(month, eta, rho) + diag_matrix(rep_vector(0.0001, K));
  L = cholesky_decompose(Sigma);
  
  // non-centered parameterization for beta
  beta = mu_beta + L * z2;
  
  // non-centered parameterization for alpha
  alpha_year = mu_alpha + sigma_alpha * z1;
  
  // linear predictor
  lambda = alpha_year + X * beta + log(offset);
}

model {
  z1 ~ normal(0, 1);
  z2 ~ normal(0, 1);
  
  rho ~ normal(0, 5);
  eta ~ normal(0, 1);

  mu_alpha ~ normal(0, 5);
  sigma_alpha ~ normal(0, 5);
  
  mu_beta ~ normal(0, 5);

  y ~ poisson_log(lambda);
}

generated quantities {
  vector[N] log_lik;
  
  for (n in 1:N) {
    log_lik[n] = poisson_log_lpmf(y[n] | lambda[n]);
  }
}
