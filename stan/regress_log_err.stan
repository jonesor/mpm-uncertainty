
data {
  int<lower=0> N;
  int<lower=0> P;
  vector[N] x;
  vector[P] xpred;
  int y[N];
  vector[N] offset;
}

parameters {
  real mu_alpha;
  real<lower=0> sigma_alpha;
  real beta;
  vector[N] z;
}

transformed parameters {
  vector[N] lambda;
  vector[N] alpha_year;

  alpha_year = mu_alpha + sigma_alpha * z;
  lambda = alpha_year + beta * x + log(offset);
}

model {
  z ~ normal(0, 1);
  
  mu_alpha ~ normal(0, 5);
  sigma_alpha ~ normal(0, 5);
  
  beta ~ normal(0, 5);
  
  y ~ poisson_log(lambda);
}

generated quantities {
  vector[N] log_lik;
  vector[N] pred_year;
  vector[P] pred;
  
  pred_year = exp(lambda - log(offset));
  pred = exp(mu_alpha + beta * xpred);
  
  for(i in 1:N) {
    log_lik[i] = poisson_log_lpmf(y[i] | lambda[i]);
  }
}
