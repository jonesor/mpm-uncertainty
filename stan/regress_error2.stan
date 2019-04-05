
data {
  int<lower=0> N;
  vector[N] x_mean;
  vector[N] x_se;
  vector[N] y_mean;
  vector[N] y_se;
}

parameters {
  real mu_alpha;
  real mu_beta;
  
  vector[N] y;
  vector[N] x;
  
  real<lower=0> sigma_y;
}

transformed parameters {
  vector[N] y_hat;
  
  for (i in 1:N) {
    y_hat[i] = mu_alpha + (mu_beta * x[i]);
  }
}

model {
  mu_alpha ~ normal(0, 1);
  mu_beta ~ normal(0, 2);
  
  sigma_y ~ normal(0, 2);
  
  y_mean ~ normal(y, y_se);
  x_mean ~ normal(x, x_se);
  
  y ~ normal(y_hat, sigma_y);
}
