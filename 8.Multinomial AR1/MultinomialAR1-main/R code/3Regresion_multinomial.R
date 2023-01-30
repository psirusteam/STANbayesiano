## Modelo binomial
library(dplyr)
library(cmdstanr)
library(lme4)

rm(list = ls())
D <- 800        # Número de dominios
K <- 3           # Número de regresoras
P <- 3           # Número de categorías
# Tiempo <- 10      # Periodos  
# phi <- -0.9       # correlación en el tiempo 

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

C <- matrix(c(1,0.7,0.7,1),2,2)
L <- chol(C)
sigma <- diag(c(1, 1))
Lambda <- sigma %*% t(L)
Z <- rbind(rnorm(D), rnorm(D))
u <- t(Lambda %*% Z)

cor(u)
sd(u[, 1])
sd(u[, 2])
mean(u[, 1])
mean(u[, 2])

### Calculando el denominador
lp1 <- crossprod(t(Xt), (beta[1, ])) + u[, 1]
lp2 <- crossprod(t(Xt), (beta[2, ])) + u[, 2]

Denominador <- 1 + exp(lp1) + exp(lp2)

# Calculando la matrix y
theta <- cbind(1 / Denominador,
               exp(lp1) / Denominador,
               exp(lp2) / Denominador)

##########################################
# Tamaño de muestra por dominio
Nd <- sample(x = 200:300, size = D, replace = TRUE ) 


y <- matrix(NA, D, P)
for (ii in 1:D) {
  y[ii, ] <- t(rmultinom(n = 1, size = Nd[ii], prob = theta[ii, ]))
}
rowSums(y) - Nd

### Validando datos simulados

sample_data <- list(D = D,
                    P = P,
                    K = K,
                    y = y,
                    X = Xt)
fit <-
  cmdstan_model(stan_file = "Modelos/Multinomial_simple.stan",
                compile = TRUE)


fit_mcmc <- fit$sample(
  num_samples = 200, 
  num_warmup = 500,
  data = sample_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4
)

fit_mcmc$print("beta")
 beta




