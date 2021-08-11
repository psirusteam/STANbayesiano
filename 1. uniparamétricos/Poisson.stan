data {
  int<lower=0> n;
  int<lower=0> y[n];
}
parameters {
  real<lower=0> theta;
}
model {
  y ~ poisson(theta);
  theta ~ gamma(38, 9);
}
