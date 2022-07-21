data {
  int<lower=0> n; // tamaño de muestra
  int<lower=0> p; // categorías
  int y[n, p]; // matriz de datos
  vector[p] alpha;
}

parameters {
  simplex[p] theta; // vector de parámetros
}

// transformed parameters {
//   real delta;
//   delta = theta[1] - theta[2];
// }

model {
  theta ~ dirichlet(alpha);
  for(i in 1:n){
    y[i, ] ~ multinomial(theta); 
  }
}

// generated quantities {
//   int ypred[k];
//   int deltapred;
//   ypred = multinomial_rng(theta, 100);
//   deltapred = ypred[1] - ypred[2];
// }
