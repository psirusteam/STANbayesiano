data {
  int<lower=1> d; // número de dominios 
  int<lower=1> p; // categorías
  int y[d, p]; // matriz de datos
  vector[p] alpha; // hiperparámetros de la prior
}

parameters {
  simplex[p] theta[d];// vector de parámetros 
}

// transformed parameters {
//   real delta;
//   delta = theta[1] - theta[2];
// }

model {
  for(i in 1:d){
    for(j in 1:p){
      target += multinomial_lpmf(y[i, ] | theta[i]); 
      target += dirichlet_lpdf(theta[j] | alpha);
    }
  }
}
