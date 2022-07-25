data {
  int<lower=1> D; // número de dominios 
  int<lower=1> P; // categorías
  int<lower=1> K; // cantidad de regresores
  int y[D, P]; // matriz de datos
  matrix[D, K] X; // matriz de covariables
}
  

parameters {
  matrix[P-1, K] beta;// matriz de parámetros 
  real<lower=0> sigma2_u;
  vector[D] u1;
}

transformed parameters {
  simplex[P] theta[D];// vector de parámetros;
  real num[D, P];
  real den[D];
  real<lower=0> sigma_u;
  
  sigma_u = sqrt(sigma2_u);
  
  for(d in 1:D){
    num[d, 1] = 1;
    for(p in 2:P){
      num[d, p] = exp(X[d, ] * beta[p-1, ]' + u1[d]) ;
    }
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
  
  u1 ~ normal(0, sigma_u);
  sigma2_u ~ inv_gamma(0.0001, 0.0001);
  
  for(p in 2:P){
    for(k in 1:K){
      beta[p-1, k] ~ normal(0, 100);
    }
  }
  
  for(d in 1:D){
    target += multinomial_lpmf(y[d, ] | theta[d, ]); 
  }

}
