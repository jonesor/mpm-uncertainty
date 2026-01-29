
data {
  int<lower=0> N;
  int<lower=0> P;
  vector[N] x;
  vector[P] xpred;
  vector[N] y;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] y_hat;
  y_hat = alpha + beta * x;
}

model {
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  sigma ~ normal(0, 10);
  
  y ~ normal(y_hat, sigma);
}

generated quantities {
  vector[N] log_lik;
  vector[N] pred_year;
  vector[P] pred;
  
  pred_year = exp(y_hat);
  pred = exp(alpha + beta * xpred);
  
  for(i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | y_hat[i], sigma);
  }
}

