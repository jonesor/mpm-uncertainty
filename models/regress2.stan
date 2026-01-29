data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}

parameters {
  real mu_alpha;
  real mu_beta;
  real<lower=0> sigma_y;
}

transformed parameters {
  vector[N] y_hat;
  
  for (i in 1:N) {
    y_hat[i] = mu_alpha + (mu_beta * x[i]);
  }
}

model {
  mu_alpha ~ normal(0, 10);
  mu_beta ~ normal(0, 10);
  sigma_y ~ normal(0, 10);
  
  y ~ normal(y_hat, sigma_y);
}
