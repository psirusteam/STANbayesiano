// The input data is a vector 'y' of length 'n'.
data {
  int<lower=0> n;             //number of experiments
  int<lower=0> m[n]; //number of observations per experiment
  int<lower=0> y[n]; //Successes per experiment
}
// The parameters accepted by the model
parameters {
  real<lower=0, upper=1> theta; //success probability
}
// The model to be estimated. We model the output
// 'y' to follow a Bernoulli distribution with mean 'theta'
model {
  // likelihood
  for(i in 1:n) {
  y ~ binomial(m, theta);
  }
  // prior
  theta ~ uniform(0, 1);  // prior for theta
}
// Posterior Predictive checks
generated quantities{
  int y_test[n];
  for(i in 1:n) {
    y_test[i] = binomial_rng(m[i],theta);
  }
}
