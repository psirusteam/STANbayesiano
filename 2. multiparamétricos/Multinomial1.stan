data {
  int<lower=0> n; // tamaño de muestra
  int<lower=0> p; // categorías
  int y[n, p]; // matriz de datos
  vector[p] alpha;
  // predicción
   int <lower = 1> ns;
   int   nd;
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

generated quantities {
  int<lower=0> ypred[ns, p]; 
  for (ii in 1:ns) {
    ypred[ii] = multinomial_rng(theta, 500);
  }
 
}
