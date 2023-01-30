### Simulación de la serie FH
library(dplyr)
library(cmdstanr)

rm(list = ls())
T_periodo <- 11 # Periodos en el tiempo
D <- 20        # Número de dominios
K <- 5         # Número de regresoras
P <- 3        # Número de categorías 

### Descomposición de Cholesky para simulación de la
### correlación de los efectos aleatorios
COV_INTRA <- diag(1,P-1)

COV_INTRA <- 0.65^abs(col(COV_INTRA)-row(COV_INTRA))   
L <- chol(COV_INTRA)
COV_Diag <- diag(1,P-1)
Lambda <- COV_Diag %*% t(L)
Z <- matrix(rnorm(D*(P-1)),nrow = P-1,ncol = D)
u <- t(Lambda %*% Z)

##################################
# Correlación temporal.
A <-  lapply(1:T_periodo,function(t){
  set.seed(1)
matrix(round(runif(n = (P-1)*(P-1),
                        min = 0.3,max = 0.5),3),
            nrow = P-1)})

##################################
## Error para el tiempo t en el área d 
sigma_Et <- diag(1,P-1)
sigma_Et <- 0.9^abs(col(sigma_Et)-row(sigma_Et))   

##################################
## ud en el tiempo t  
U <- lapply(1:T_periodo, function(x){
  matrix(rnorm(D * (P-1), 5, 1.5),
              D, (P-1))
})

##################################
# Coeficiente del modelo
beta1 <- c(0.2, -.4)
beta <- matrix(beta1, P-1, K,byrow = TRUE) # coefientes del modelo.
##################################
# Covariables
Xt <- cbind(1, matrix(rnorm(D * (K - 1), 5, 3),
                      D, (K - 1)))

##################################
# Correlación temporal. 
Udt <- lapply(2:T_periodo, function(tt){
crossprod(t(U[[tt-1]]),A[[tt]]) +  
    mvtnorm::rmvnorm(n = D, sigma = sigma_Et) 
})

### Calculando el denominador
lp <- lapply(1:(T_periodo - 1), function(tt) {
  lapply(2:P, function(pp){
    crossprod(t(Xt), (beta[pp-1, ])) 
  })%>% Reduce(f = cbind, x = .) + Udt[[tt]] +u
})


theta <- lapply(1:(T_periodo - 1), function(tt) {
  lptemp <- lp[[tt]]
  paso <- matrix(NA, D, P)
  Denominador <- 1 + rowSums(exp(lptemp))
  paso[, 1] <- 1 / Denominador
  for (pp in 2:P) {
    paso[, pp] <- exp(lptemp[, pp - 1]) / Denominador
  }
  round(paso,4)
})
lapply(theta, rowSums)


##################################
# Datos 
sample_data <- data.frame(id= 1:D, Xt)

Y <- array(NA, dim = c(D,P,T_periodo-1))

for(tt in 1:(T_periodo-1)){
   #variables regresoras
  Nd <- sample(x = 200:300,size = D,replace = TRUE)
  yd <- matrix(NA, D, P)
  colnames(yd) <- paste0("y", 1:P)
    for (dd in 1:D) {
    yd[dd, ] <- t(rmultinom(n = 1, size = Nd[dd],
                           prob = theta[[tt]][dd, ]))
  }
  
  Y[,,tt] <- yd
}

array_Udt <- array(NA, dim = c(D,P-1,T_periodo-1))
for(tt in 1:(T_periodo-1)){
  array_Udt[,,tt] <- Udt[[tt]]
}

sample_data2 <- list(
  D = D,                       ## Número de dominios
  P = K,                       ## Número de regresoras
  K = P,                       ## Número de categoría
  T = T_periodo - 1,           ## Número de regresoras
  Yd = Y,                      ## Vector con el número de éxitos.
                               ## Vector con el número de ensayos.
  Nd = t(apply(Y, MARGIN = 1, colSums)) %>%
    array(data = ., dim = c(D, 1, T_periodo - 1)),
  
  Udt = array_Udt,             ## Vector con el número de éxitos.
  X = (Xt)                     ## Matriz de regresores
)


file.remove("Modelos/Regresion_temp.exe")
fit <-
  cmdstan_model(stan_file = "Modelos/Regresion_temp.stan",
                compile = TRUE)

fit_mcmc <- fit$sample(
  num_samples = 10,
  num_warmup  = 10,
  data = sample_data2,
  seed = 123,
  chains = 1,
  parallel_chains = 1
)

betaE <- fit_mcmc$summary("beta") %>% select(mean) 
matrix(data = betaE$mean,ncol = K)
beta
fit_mcmc$summary("rho")
phi
sigma2e <- fit_mcmc$summary("sigma2_e") %>% select(mean)
sqrt(sigma2e$mean)[-1]
sigma_Et[-(1:2)]
sigmav <- fit_mcmc$summary("sigma_v")%>% arrange(desc(sd)) %>% 
  select(variable:sd)
sigmav %>% data.frame()
sigmav$mean
sigma_v
ypred <- fit_mcmc$summary("Ydpred")
matrix(ypred$mean,nrow = D) -Y

library(bayesplot)
library(patchwork)
(mcmc_dens_chains(fit_mcmc$draws("rho")) +
    mcmc_areas(fit_mcmc$draws("rho")))/ 
  mcmc_trace(fit_mcmc$draws("rho"))

mcmc_dens_chains(fit_mcmc$draws(c("sigma2_e[2]","sigma2_e[3]","sigma2_e[4]","sigma2_e[5]","sigma2_e[6]")))
mcmc_areas(fit_mcmc$draws(c("sigma2_e[2]","sigma2_e[3]","sigma2_e[4]","sigma2_e[5]","sigma2_e[6]"))) 
mcmc_trace(fit_mcmc$draws(c("sigma2_e[2]","sigma2_e[3]","sigma2_e[4]","sigma2_e[5]","sigma2_e[6]")))

mcmc_dens_chains(fit_mcmc$draws("sigma_v"))
mcmc_areas(fit_mcmc$draws("sigma_v")) 
mcmc_trace(fit_mcmc$draws("sigma_v"))

mcmc_dens_chains(fit_mcmc$draws("beta"))
mcmc_areas(fit_mcmc$draws("beta")) 
mcmc_trace(fit_mcmc$draws("beta"))
