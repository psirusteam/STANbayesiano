// Modelo multinomial con ud fijo para todo tiempo t
// El resultado es una matrix de bestas de K*(P-1)

data {
  int<lower=1> D; // número de dominios 
  int<lower=1> P; // categorías
  int<lower=1> K; // cantidad de regresores
  int<lower=2> T; // Tiempo
  int y[D, P, T]; // matriz de datos
  matrix[D, K] X; // matriz de covariables
}



parameters {
 matrix[K,P-1] beta;// matriz de parámetros 
   vector<lower=0>[P-1] sigma_u;       // random effects standard deviations
  // declare L_u to be the Choleski factor of a 2x2 correlation matrix
  cholesky_factor_corr[P-1] L_u;
  matrix[P-1, D] z_u;  
}

transformed parameters {
  matrix[D,1]lp1[T];
  matrix[D,1]lp2[T];
  matrix[D,1]Den[T];
  simplex[P] theta[D,T];
  matrix[P-1, D] u; // random effect matrix
  u = diag_pre_multiply(sigma_u, L_u) * z_u;
  
  
  for(tt in 1:T){
    lp1[tt,,1] = X*beta[,1] +  to_vector(u[1, ]); 
    lp2[tt,,1] = X*beta[,2] +  to_vector(u[2, ]); 
    Den[tt,,1] = 1 + exp( lp1[tt,,1]) + exp( lp2[tt,,1]);
   
    for(dd in 1:D){
    theta[dd,tt,1] =  1/Den[tt,dd,1];
    theta[dd,tt,2] =  exp( lp1[tt,dd,1])/Den[tt,dd,1];
    theta[dd,tt,3] =   exp( lp2[tt,dd,1])/Den[tt,dd,1];
    }

}
  
}
model {
  to_vector(beta) ~ normal(0,1000);
  L_u ~ lkj_corr_cholesky(1); // LKJ prior for the correlation matrix
  to_vector(z_u) ~ normal(0, 1);
  sigma_u ~ cauchy(0, 50);
  
  for (tt in 1:T){
    for(dd in 1:D){
    target += multinomial_lpmf(y[dd,, tt] | theta[dd,tt, ]); 
  }
  }
}

generated quantities {
  matrix[2, 2] Omega;
  vector<lower=0>[2] sdcomprobar;
  sdcomprobar[1] = sd(u[1, ]);
  sdcomprobar[2] = sd(u[2, ]);

  Omega = L_u * L_u'; // so that it return the correlation matrix
}

