data {
  int<lower=1> T; // Tiempo
  int<lower=1> K; // Regresoras 
  int<lower=0> D; // Dominios 
  matrix[D, K] X; // Matriz [D,P]
  int <lower=0> Yd[D,T]; // array [D,1,T]
  int <lower=0> Nd[D]; // array [D,1,T]
  matrix [D,T] udt; // Matriz [D,T].
}


parameters {
  matrix [K,T] beta;
  real <lower=0> sigma2_ud;
  real <lower=0> sigma2_ut[T];
  vector[D]ud;
  real<lower=-1,upper=1> phi;
  
}

transformed parameters {
  matrix [D,T] lp ;
  matrix [D,T] eta ;
  real <lower=0> sigma_ud;
  real <lower=0> sigma_ut[T];
  sigma_ud = sqrt(sigma2_ud);
  sigma_ut = sqrt(sigma2_ut);
  
  for(tt in 1:T ){
    lp[,tt] = X*beta[,tt] + ud +  sigma_ut[tt]*udt[,tt];
    eta[,tt] = inv_logit(lp[,tt]);
  }
  
}

model {
  to_vector(beta) ~ normal(0,1000);
  sigma2_ud ~ inv_gamma(0.001, 0.001);
  sigma2_ut ~ inv_gamma(0.001, 0.001);
  ud ~ normal(0, sigma_ud);
  phi ~ uniform(-1, 1);
  
  for(tt in 2:T){
  udt[,tt] ~ normal(phi*udt[,tt-1], sqrt(1 - pow(phi, 2) )); 
  }

  
  
  for(tt in 1:T){
    Yd[,tt] ~ binomial(Nd, eta[,tt]); 
  }
}

