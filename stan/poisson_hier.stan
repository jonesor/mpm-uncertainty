
data {
  int<lower=0> J;           // # of groups
  int<lower=0> J_unobs;     // # of unobserved groups to simulate
  int<lower=0> y[J];        // count
  vector[J] offset;         // offset
}

parameters {
  vector[J] z;
  real mu;
  real<lower=0> sigma;
}

transformed parameters {
  vector[J] alpha;
  vector[J] lambda;
  
  alpha = mu + z * sigma;
  lambda = alpha + offset;
}

model {
  z ~ normal(0, 1);
  mu ~ normal(0, 2);
  sigma ~ normal(0, 1);
  
  y ~ poisson_log(lambda);
}

generated quantities {
  vector[J_unobs] alpha_new;
  
  for(j in 1:J_unobs) {
    alpha_new[j] = normal_rng(mu, sigma);
  }
}
