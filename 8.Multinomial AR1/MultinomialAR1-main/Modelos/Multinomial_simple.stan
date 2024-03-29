data {
  int<lower=1> D; // número de dominios 
  int<lower=1> P; // categorías
  int<lower=1> K; // cantidad de regresores
  int y[D, P]; // matriz de datos
  matrix[D, K] X; // matriz de covariables
}
  

parameters {
  matrix[P-1, K] beta;// matriz de parámetros 
  vector<lower=0>[P-1] sigma_u;       // random effects standard deviations
  // declare L_u to be the Choleski factor of a 2x2 correlation matrix
  cholesky_factor_corr[P-1] L_u;
  matrix[P-1, D] z_u;                  
}

transformed parameters {
  simplex[P] theta[D];// vector de parámetros;
  real num[D, P];
  real den[D];
  // this transform random effects so that they have the correlation
  // matrix specified by the correlation matrix above
  matrix[P-1, D] u; // random effect matrix
  u = diag_pre_multiply(sigma_u, L_u) * z_u;
  
  for(d in 1:D){
    num[d, 1] = 1;
    num[d, 2] = exp(X[d, ] * beta[1, ]' + u[1, d]) ;
    num[d, 3] = exp(X[d, ] * beta[2, ]' + u[2, d]) ;
    
    den[d] = sum(num[d, ]);
  }
  
  for(d in 1:D){
    for(p in 2:P){
    theta[d, p] = num[d, p]/den[d];
    }
    theta[d, 1] = 1/den[d];
  }
}

model {
  L_u ~ lkj_corr_cholesky(1); // LKJ prior for the correlation matrix
  to_vector(z_u) ~ normal(0, 1);
  sigma_u ~ cauchy(0, 50);
  
  for(p in 2:P){
    for(k in 1:K){
      beta[p-1, k] ~ normal(0, 100);
    }
    }
  
  for(d in 1:D){
    target += multinomial_lpmf(y[d, ] | theta[d, ]); 
  }
}

  
generated quantities {
  matrix[2, 2] Omega;
  vector<lower=0>[2] sdcomprobar;
  sdcomprobar[1] = sd(u[1, ]);
  sdcomprobar[2] = sd(u[2, ]);

  Omega = L_u * L_u'; // so that it return the correlation matrix
}
