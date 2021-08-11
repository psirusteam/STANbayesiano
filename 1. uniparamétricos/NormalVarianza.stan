data {
  int<lower=0> n;
  real mu;
  real y[n];
}
parameters {
  real<lower = 0> tau;
}
transformed parameters {
  real tau2;
  tau2 = pow(tau, 2);
}
model {
  y ~ normal(mu, tau);
  //tau ~ inv_gamma(0.1, 0.1);
  //tau ~ uniform(1, 30);
  //tau ~ inv_gamma(1, 1);
  tau ~ cauchy(0, 25);
}
generated quantities {
  real y_test[n];
  for(i in 1:n) {
    y_test[i] = normal_rng(mu, tau);
  }
}
