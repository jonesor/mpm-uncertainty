data {
  int<lower=1> T;
  array[T] int<lower=0> recruits;
  vector[T] x;
  vector<lower=0>[T] n_repro;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma_proc;
  vector[T] eta_raw;
}

transformed parameters {
  vector[T] eta;

  eta = alpha + beta * x + sigma_proc * eta_raw;
}

model {
  alpha ~ normal(0, 2);
  beta ~ normal(0, 1);
  sigma_proc ~ exponential(1);
  eta_raw ~ std_normal();
  recruits ~ poisson_log(log(n_repro) + eta);
}

generated quantities {
  vector[T] log_lik;

  for (t in 1:T) {
    log_lik[t] = poisson_log_lpmf(recruits[t] | log(n_repro[t]) + eta[t]);
  }
}
