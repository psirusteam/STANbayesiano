data {
  int<lower=1> n;                    // sample size
  int<lower=1> N;                    // population size
  int<lower=1> p;                     // p predictors
  matrix[n, p] Xsam;         // design matrix
  matrix[N, p] Xpred;        // predictive matrix
  vector<lower=0>[n] y;      // response variable
}

parameters {
  vector[p] beta; //the regression parameters
  real<lower=0> phi; //the variance parameter
}

transformed parameters {
  vector[n] LPsam; //linear predictor
  vector[n] theta; //the expected values 
  real<lower=0> a; // shape parameter for the gamma distribution
  vector[n] b; //rate parameter for the gamma distribution

  a = 1 / phi; // shape parameter for the gamma distribution
  LPsam = Xsam * beta; // linear predictor
  for (i in 1:n) { 
    theta[i] = exp(LPsam[i]); 
    b[i] = a / (theta[i]);
  }
    
}

model {
 // priors
  beta ~ normal(0, 100);
  phi ~ inv_gamma(0.0001, 0.0001);
  // likelihood
  y ~ gamma(a, b);
}

generated quantities {
  vector[n] ysam;
  vector[N] LPpred;
  vector[N] ypred;
  vector[N] thetapred;
  vector[N] bpred; 
  
  for (i in 1:n) {
    ysam[i] = gamma_rng(a, b[i]); //posterior draws to get posterior predictive checks
  }

  LPpred = Xpred * beta;
  for (j in 1:N) {
    thetapred[j] = exp(LPpred[j]);
    bpred[j] = a / thetapred[j];
    ypred[j] = gamma_rng(a, bpred[j]);
  }
}

  
