data {
  int<lower=0> n;
  vector[n] y;
}
parameters {
  real<lower=0> theta;
}
transformed parameters {
  real<lower=0> invtheta = 1/theta;
}
model {
  y ~ exponential(theta);
  theta ~ gamma(0.1, 0.1);
}
