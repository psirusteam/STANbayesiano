data {
  int<lower=1> D; // número de dominios 
  int<lower=1> P; // categorías
  int<lower=1> K; // cantidad de regresores
  int y[D, P]; // matriz de datos
  matrix[D, K] X; // matriz de covariables
}
  

parameters {
  matrix[P-1, K] beta;// matriz de parámetros 
  real<lower=0> sigma2_u1;
  real<lower=0> sigma2_u2;
  vector[D] u1;
  vector[D] u2;
}

transformed parameters {
  simplex[P] theta[D];// vector de parámetros;
  real num[D, P];
  real den[D];
  real<lower=0> sigma_u1;
  real<lower=0> sigma_u2;
  
  sigma_u1 = sqrt(sigma2_u1);
  sigma_u2 = sqrt(sigma2_u2);
  
  for(d in 1:D){
    num[d, 1] = 1;
    num[d, 2] = exp(X[d, ] * beta[1, ]' + u1[d]) ;
    num[d, 3] = exp(X[d, ] * beta[2, ]' + u2[d]) ;
    
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
  
  u1 ~ normal(0, sigma_u1);
  //sigma2_u1 ~ inv_gamma(0.0001, 0.0001);
  sigma2_u1 ~  cauchy(0, 50);
   
  u2 ~ normal(0, sigma_u2);
  //sigma2_u2 ~ inv_gamma(0.0001, 0.0001);
  sigma2_u2 ~  cauchy(0, 50);
  
  for(p in 2:P){
    for(k in 1:K){
      beta[p-1, k] ~ normal(0, 100);
    }
  }
  
  for(d in 1:D){
    target += multinomial_lpmf(y[d, ] | theta[d, ]); 
  }

}
