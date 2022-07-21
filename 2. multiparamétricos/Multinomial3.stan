data {
  int<lower=0> d; // número de dominios 
  int<lower=0> p; // categorías
  int y[d, p]; // matriz de datos
  vector[p] alpha;
}

parameters {
  matrix [d,p] theta;// vector de parámetros 

}

// transformed parameters {
//   real delta;
//   delta = theta[1] - theta[2];
// }

model {
  for(i in 1:d){
   theta[i,] ~ dirichlet(alpha);
   y[i, ] ~ multinomial(theta[i,]); 
  }
}

// generated quantities {
//   int ypred[k];
//   int deltapred;
//   ypred = multinomial_rng(theta, 100);
//   deltapred = ypred[1] - ypred[2];
// }
