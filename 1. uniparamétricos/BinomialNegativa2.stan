data {
  int<lower=0> n;
  int<lower=0> k[n];
  int<lower=0> y[n];
}
transformed data {
  int<lower=0> m[n];
  for(i in 1:n){
    m[i] = y[i] - k[i];
  }
}
parameters {
  real<lower=0> b;
}
transformed parameters {
  real<lower=0> theta;
  theta = b/(b + 1);
}
model {
  for(i in 1:n){
    m[i] ~ neg_binomial(k[i], b);
  }
  theta ~ beta(0.5, 0.5);
}

