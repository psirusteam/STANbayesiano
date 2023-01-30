## Modelo binomial
library(dplyr)
library(cmdstanr)
library(lme4)

rm(list = ls())
D <- 400        # Número de dominios
K <- 3           # Número de regresoras
Tiempo <- 10      # Periodos  
phi <- -0.9       # correlación en el tiempo 

## 

##################################
# Covariables
Xt <- cbind(runif(D), runif(D),  runif(D))
###################################
betas <- matrix(c(1, -2, 3), 3)

###################################
# Error área 
ud <- rnorm(D, 0,sd = 0.5)

###################################
# Error en el tiempo
ut <- matrix(NA,  D, Tiempo)
ut[,1] <- rnorm(D, 0,sd = sqrt(1-phi^2))

for(tt in 2:(Tiempo)){
  et <- rnorm(D, 0, sd = sqrt(1-phi^2))
  ut[,tt] <- phi*ut[,tt-1] +  et 
}

###################################
# Parte lineal. 
lp <- matrix(NA, nrow = D, ncol = Tiempo)
eta <- matrix(NA, nrow = D, ncol = Tiempo)

for(tt in 1:Tiempo){
  lp[,tt] <- Xt %*% betas  + ud + ut[,tt]
  eta[,tt] <- boot::inv.logit(lp[,tt])  
}
summary(eta)
##################################
# Número de ensayos 
Nd <- sample(x = 200:300,size = D, replace = TRUE)

##################################
# Creando yd Numero de éxitos. 
yd <- matrix(NA, nrow = D, ncol = Tiempo)

for( tt in 1:Tiempo){
  yd[,tt] <- rbinom(n = D, size = Nd, prob = eta[,tt])  
}

# summary(yd/Nd - eta)
##################################
# Creando yd Numero de éxitos. 
id_Dominio <- 1:D


## Validando el modelo 

for(tt in 1:Tiempo){
  fit <- glmer(cbind(yd[,tt], Nd-yd[,tt]) ~ 0 + (1|id_Dominio) +Xt,
               family = binomial(link = "logit"))
  print(fit)
}

pacf(ts(as.numeric(yd[3,])))


sample_data2 <- list(
  D = D,           ## Número de dominios
  K = K,           ## Número de regresoras
  T = Tiempo, ## Número de regresoras
  Yd =  yd,        ## Vector con el número de éxitos.
  Nd = Nd,         ## Vector con el número de ensayos.
  udt = (ut),      ## Vector con el número de éxitos.
  X = (Xt)         ## Matriz de regresoras
)


#file.remove("Modelos/Regresion_logit_binomial_AR1_V0.exe")
fit <-
  cmdstan_model(stan_file = "Modelos/Regresion_logit_binomial_AR1.stan",
                compile = TRUE)

fit_mcmc <- fit$sample(
  num_samples = 100,
  num_warmup  = 100,
  data = sample_data2,
  seed = 123,
  chains = 4,
  parallel_chains = 4
)

betaE <- fit_mcmc$summary("beta") %>% select(mean) 
matrix(data = betaE$mean, ncol = K, byrow = TRUE)
betas
fit_mcmc$summary("sigma_ud")

fit_mcmc$summary("phi")
phi

fit_mcmc$summary("sigma_ut")
