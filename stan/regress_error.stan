
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y_mu;
  vector[N] y_se;
}

parameters {
  vector[N] z;
  real alpha;
  real beta;
  real<lower=0> tau;
}

transformed parameters {
  vector[N] theta;
  vector[N] y_hat;

  y_hat = alpha + beta * x;
  theta = y_hat + tau * z;
}

model {
  z ~ normal(0, 1);
  tau ~ normal(0, 10);
  
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  
  y_mu ~ normal(theta, y_se);
}
