## Modelo binomial
library(dplyr)
library(cmdstanr)
library(lme4)

rm(list = ls())
D <- 200        # Número de dominios
K <- 3           # Número de regresoras
P <- 3           # Número de categorías
Tiempo <- 10      # Periodos  
phi <- 0.9       # correlación en el tiempo 
A <- diag(phi, P-1)
## 

##################################
# Covariables
Xt <- cbind(runif(D), runif(D),  runif(D))
###################################
# coeficientes del modelo.
beta <- matrix(c(-2, -2, 1, 1), P-1, K) 

##################################
## Error del área 
##################################

### Descomposición de Cholesky para simulación de la
### correlación de los efectos aleatorios
ut <- array(NA, dim = c(D, P-1,Tiempo))

C <- matrix(c(1,0.7,0.7,1),2,2)
L <- chol(C)
sigma <- diag(c(1, 1))
Lambda <- sigma %*% t(L)

Z <- rbind(rnorm(D), rnorm(D))
ut[,,1] <- t(Lambda %*% Z)

Sigma_e = diag(1,P-1) - A %*% t(A)  

for(tt in 2:Tiempo){
  et <- mvtnorm::rmvnorm(n = D, mean = c(0,0), sigma = Sigma_e)
    ut[,,tt] <- (ut[,,tt-1])%*%A  + et  
}


### Calculando el denominador
lpt1 <- array(NA, dim = c(D, 1,Tiempo))
lpt2 <- array(NA, dim = c(D, 1,Tiempo))
Denominador <- array(NA, dim = c(D, 1,Tiempo))
theta <- array(NA, dim = c(D, P,Tiempo))

for(tt in 1:Tiempo){
lpt1[,,tt] <- crossprod(t(Xt), (beta[1, ])) + ut[, 1,tt]
lpt2[,,tt] <- crossprod(t(Xt), (beta[2, ])) + ut[, 2,tt]
Denominador[,,tt] <- 1 + exp(lpt1[,,tt]) + exp(lpt2[,,tt])

# Calculando la matrix de theta
theta[,,tt] <- cbind(1 / Denominador[,,tt],
               exp(lpt1[,,tt]) / Denominador[,,tt],
               exp(lpt2[,,tt]) / Denominador[,,tt])

}


##########################################
# Tamaño de muestra por dominio
Nd <- sample(x = 200:300, size = D, replace = TRUE ) 

y <- array(NA, dim = c(D, P,Tiempo))
for(tt in 1:Tiempo){
for (ii in 1:D) {
  y[ii, ,tt ] <- t(rmultinom(n = 1, size = Nd[ii], prob = theta[ii,, tt]))
}
}  


### Validando datos simulados

sample_data <- list(D = D,
                    P = P,
                    K = K,
                    T = Tiempo,
                    y = y,
                    X = Xt,
                    ut = ut)
fit <-
  cmdstan_model(stan_file = "Modelos/Regresion_logit_multnomial_AR1_rho_V02.stan",
                compile = TRUE)


fit_mcmc <- fit$sample(
  num_samples = 1000, 
  num_warmup =  1000,
  data = sample_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 20,
)

fit_mcmc$print("beta")
beta

fit_mcmc$print("sigma_u")

fit_mcmc$print("Omega")
fit_mcmc$print("A") 
A

