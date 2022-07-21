
data {
  int<lower=0> k;
  int y[k];
  vector[k] alpha;
}

parameters {
  simplex[k] theta;
}

transformed parameters {
  real delta;
  delta = theta[1] - theta[2];
}

model {
  y ~ multinomial(theta);
  theta ~ dirichlet(alpha);
}

// generated quantities {
//   int ypred[k];
//   int deltapred;
//   ypred = multinomial_rng(theta, 100);
//   deltapred = ypred[1] - ypred[2];
// }
