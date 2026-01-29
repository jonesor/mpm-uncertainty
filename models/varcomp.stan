data {
  int N;
  vector[N] y_mean;
  vector[N] y_se;
}

transformed data {
  vector[N] y_var;
  real var_w;
  
  for (i in 1:N) y_var[i] = pow(y_se[i], 2.0);
  var_w = mean(y_var);
}

parameters {
  real mu;                // population treatment effect
  real<lower=0> tau;      // standard deviation in treatment effects
  vector[N] eta;          // unscaled deviation from mu by school
}
transformed parameters {
  vector[N] theta = mu + tau * eta;        // school treatment effects
}
model {
  eta ~ normal(0, 1);       // prior log-density
  y_mean ~ normal(theta, y_se); // log-likelihood
}

generated quantities {
  real var_a = pow(tau, 2.0);
  real pvar_w = var_w / (var_w + var_a);
}
