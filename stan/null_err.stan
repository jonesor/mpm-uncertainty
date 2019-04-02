
data {
  int<lower=0> N;
  int y[N];
  vector[N] offset;
}

parameters {
  real mu_alpha;
  real<lower=0> sigma_alpha;
  vector[N] z;
}

transformed parameters {
  vector[N] alpha_year;
  vector[N] lambda;
  
  // non-centered parameterization for alpha
  alpha_year = mu_alpha + sigma_alpha * z;
  
  // linear predictor
  lambda = alpha_year + log(offset);
}

model {
  z ~ normal(0, 1);
  
  mu_alpha ~ normal(0, 5);
  sigma_alpha ~ normal(0, 5);
  
  y ~ poisson_log(lambda);
}

generated quantities {
  vector[N] log_lik;

  for (n in 1:N) {
    log_lik[n] = poisson_log_lpmf(y[n] | lambda[n]);
  }
}
