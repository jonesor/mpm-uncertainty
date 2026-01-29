
data { 
  int<lower=0> J;             // # of groups
  int<lower=0> J_unobs;       // # of unobserved groups to simulate
  int<lower=0> K;             // # of possible outcomes
  int<lower=0> N;             // # total observations
  int<lower=0> group[N];      // group id for observation n
  int<lower=1,upper=K> y[N];  // outcome for observation n
}

parameters {
  vector[K-1] mu;               // pop. mean log-odds {hit, walk, other}
  vector<lower=0>[K-1] sigma;   // pop. sd log-odds {hit, walk, other}
  matrix[K-1, J] eta;           // non-centered prior
}

transformed parameters {
  matrix[K, J] alpha;         // group-specific log-odds {hit, walk, other}
  
  for(k in 1:(K-1)) {
    alpha[k,] = mu[k] + sigma[k] * row(eta, k);
  }
  
  alpha[K,] = rep_row_vector(0.0, J);
}

model {
  mu ~ normal(0, 1);
  sigma ~ normal(0, 1);
  to_vector(eta) ~ normal(0, 1);
  
  for(n in 1:N) {
    y[n] ~ categorical_logit(col(alpha, group[n]));
  }
}

generated quantities {
  matrix[K-1, J] theta;   // transform alpha to probability scale
  matrix[K-1, J_unobs] theta_new;
  matrix[K, J_unobs] alpha_new;
  
  
  vector[K] mu_full;
  vector[K] mu_soft;
  
  mu_full = rep_vector(0.0, K);
  mu_full[1:(K-1)] = mu;
  mu_soft = softmax(mu_full);
  
  for(j in 1:J) {
    theta[,j] = softmax(col(alpha, j))[1:(K-1)];
  }
  
  for(j in 1:J_unobs) {
    for(k in 1:(K-1)) {
      alpha_new[k, j] = normal_rng(mu[k], sigma[k]);
    }
    alpha_new[K, j] = 0;
    theta_new[,j] = softmax(col(alpha_new, j))[1:(K-1)];
  }
}
