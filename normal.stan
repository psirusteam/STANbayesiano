data {
  int<lower = 1> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real a;                // intercept
  real b;                // slope
  real<lower = 0> sigma; // noise standard deviation
}
model {
  y ~ normal(a + b * x, sigma); // likelihood
}
generated quantities {
  vector[N] y_test;
  for(i in 1:N) {
    y_test[i] = normal_rng(a + b * x[i], sigma);
  }
}
