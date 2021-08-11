data {
  int<lower=0> n;
  real y[n];
}
parameters {
  real theta;
}
model {
  y ~ normal(theta, 0.1);
  theta ~ normal(2.8, 0.23);
}
