functions {
  matrix pred_theta(matrix Xp, matrix Zp, int p, matrix beta, matrix u){
  int D1 = rows(Xp);
  real num1[D1, p];
  real den1[D1];
  matrix[D1,p] theta_p;

for(d in 1:D1){ 
  for(pp in 1:p){
 theta_p[d,pp] = inv_logit(Xp[d, ] * beta[pp,]' + Zp[d, ] * u[pp,]');
 }
 }  
  
  return theta_p  ;
  }
  
}



data {
  int<lower=1> D;    // número de postestrto 
 // int<lower=1> D1;   // número de dominios por predesir 
  int<lower=1> P;    // categorías
  int<lower=1> K;  // cantidad de regresores
  int<lower=1> Kz; // cantidad de regresores en Z
  int y[D, P];       // matriz de datos
  matrix[D, K] X; // matriz de covariables
  matrix[D, Kz] Z; // matriz de covariables
 // matrix[D1, K] Xp; // matriz de covariables
 int<lower=1> nd[D];
}
  

parameters {
  matrix[P, K] beta;// matriz de parámetros 
  vector<lower=0>[P] sigma_u;       // random effects standard deviations
  // declare L_u to be the Choleski factor of a 2x2 correlation matrix
  cholesky_factor_corr[P] L_u;
  matrix[P, Kz] z_u;                  
}

transformed parameters {
  matrix [D,P] theta;// vector de parámetros;
  // this transform random effects so that they have the correlation
  // matrix specified by the correlation matrix above
  matrix[P, Kz] u; // random effect matrix
  u = diag_pre_multiply(sigma_u, L_u) * z_u;
  
for(d in 1:D){ 
  for(pp in 1:P){
 theta[d,pp] = inv_logit(X[d, ] * beta[pp,]' + Z[d, ] * u[pp,]');
 }
 }
}

model {
  L_u ~ lkj_corr_cholesky(1); // LKJ prior for the correlation matrix
  to_vector(z_u) ~ normal(0, 1);
  to_vector(sigma_u) ~ cauchy(0, 50);
  to_vector(beta) ~ normal(0, 100);
 for(d in 1:D){
   for(pp in 1:P){
  y[d,pp] ~ binomial(nd[d], theta[d,pp]);
}
}
}


generated quantities {
  // predict
  // matrix[D1,P] theta_p;// vector de parámetros;
  matrix[P, P] Omega;
  vector<lower=0>[P] sdcomprobar;
  sdcomprobar[1] = sd(u[1, ]);
  sdcomprobar[2] = sd(u[2, ]);
  sdcomprobar[3] = sd(u[3, ]);
  sdcomprobar[4] = sd(u[4, ]);
  sdcomprobar[5] = sd(u[5, ]);
  sdcomprobar[6] = sd(u[6, ]);
  sdcomprobar[7] = sd(u[7, ]);
  sdcomprobar[8] = sd(u[8, ]);

  Omega = L_u * L_u'; // so that it return the correlation matrix
// predicción

//theta_p = pred_theta(Xp,P, beta) ;

}
