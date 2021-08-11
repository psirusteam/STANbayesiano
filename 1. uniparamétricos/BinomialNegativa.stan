data {
  int<lower=0> k;
  int<lower=0> y;
}
transformed data {
  int<lower=0> m;
  m = y - k;
}
parameters {
  real<lower=0> beta;
}
transformed parameters {
  real<lower=0> theta;
  theta = beta/(beta + 1);
}
model {
  m ~ neg_binomial(k, beta);
  theta ~ beta(0.5, 0.5);
}

