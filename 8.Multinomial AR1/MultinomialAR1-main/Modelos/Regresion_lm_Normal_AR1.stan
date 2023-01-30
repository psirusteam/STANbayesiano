data {
  int<lower=1> T; // Tiempo
  int<lower=1> P; // Regresoras 
  int<lower=0> D; // Dominios 
  vector[P] X[D]; // Matriz [D,P]
  vector[T] Y[D]; // Matriz [D,T]
  vector[D]Udt[T]; // Matriz [D,T].
  cov_matrix[T] SIGMA; // Matriz [T,T]
}

parameters {
  matrix[T,P] beta;
  real<lower=-1,upper=1> rho;
  real<lower=0> sigma2_e[T];
  real<lower=0> sigma2_v[D];
  vector[D] vd ;  

}

transformed parameters {
real<lower=0> sigma_v[D];
  vector[T] mu[D];
sigma_v = sqrt(sigma2_v);

// Resultados por dominio 

for (nn in 1:D){
    mu[nn] = beta * X[nn] + vd[nn] + to_vector(Udt[,nn]);
}

}
model {
  to_vector(sigma2_e)  ~ inv_gamma(0.01,0.01);
  sigma2_v  ~ inv_gamma(0.01,0.01);
  vd ~ normal(0,sigma_v);
  to_vector(beta) ~ normal(0,1000);
 // tiempo 
for(tt in 2:T){
Udt[tt] ~ normal(rho*Udt[tt-1],sqrt(sigma2_e[tt])); 
}

  Y ~ multi_normal(mu, SIGMA); 
 }
