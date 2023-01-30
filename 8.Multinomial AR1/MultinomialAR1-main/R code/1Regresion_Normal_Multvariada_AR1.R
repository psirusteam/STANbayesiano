### Simulación de la serie FH
library(dplyr)
library(cmdstanr)
library(forecast)

rm(list = ls())
T_periodo <- 10 # Periodos en el tiempo
phi <- -0.3     # Correlación temporal.
D <- 50        # Número de dominios
K <- 3         # Número de regresoras

##################################
## Error muestral en los t tiempos
Sigmadt <- matrix(NA,T_periodo,T_periodo)
## Sigmadt es conocida. 
Sigmadt <- 0.8^abs(col(Sigmadt)-row(Sigmadt))   
edt <- mvtnorm::rmvnorm(n = D, sigma = Sigmadt)

##################################
## Error para el tiempo t en el área d 
sigma_Et <- 2
Edt <-  matrix(rnorm(D*T_periodo,0,sigma_Et),D,T_periodo)
sigma_Et <- apply(Edt, 2,sd)
##################################
## Error en el área d                
sigma_v <- runif(D)
vd <- rnorm(D,0,sigma_v) 

##################################
# Coeficiente del modelo
beta <- c(2, 1, 0.3)

##################################
# Covariables
Xt <- cbind(1, matrix(rnorm(D * (K - 1), 5, 3),
                      D, (K - 1)))
##################################
## ud en el tiempo t  
u <- matrix(NA,
            D, T_periodo)
u[,1] <- rnorm(D)
##################################
# Datos 
sample_data <- data.frame(id= 1:D, Xt)
#colnames(Xt) <- paste0("X",1:K,"t")

for (tt in 2:T_periodo) {
  #variables regresoras
  u[,tt] <- phi*u[,tt-1] + Edt[,tt]

  yt <- crossprod(t(Xt), (beta)) + vd + u[,tt] + edt[,tt]
  temp <- data.frame(yt, ut = u[,tt])
  colnames(temp) <- paste0(colnames(temp), tt - 1)
  sample_data <- bind_cols(sample_data, temp)
}

Y <- sample_data %>% select(matches("yt")) 
matplot(t(Y),type = "l")
udt <- sample_data %>% select(matches("ut")) 

acf(ts(as.numeric(Y[1,])))

auto.arima(ts(as.numeric(Y[1,])))

sample_data2 <- list(
  D = D,           ## Número de dominios
  P = K,           ## Número de regresoras
  T = T_periodo-1,           ## Número de regresoras
  SIGMA = Sigmadt[-1,-1],      ## Matriz de correlación
  Y = (Y),               ## Vector con el número de exitos.
  Udt = t(udt),          
  X = (Xt)            ## Matriz de regresores
)

file.remove("Modelos/Regresion_lm_Normal_AR1.exe")
fit <-
  cmdstan_model(stan_file = "Modelos/Regresion_lm_Normal_AR1.stan",
                compile = TRUE)

fit_mcmc <- fit$sample(
  num_samples = 500,
  num_warmup = 500,
  data = sample_data2,
  seed = 123,
  chains = 4,
  parallel_chains = 4
)

betaE <- fit_mcmc$summary("beta") %>% select(mean) 
matrix(data = betaE$mean,ncol = K)

fit_mcmc$summary("rho")
phi
sigma2e <- fit_mcmc$summary("sigma2_e") %>% select(mean)
sqrt(sigma2e$mean)[-1]
sigma_Et[-(1:2)]
sigmav <- fit_mcmc$summary("sigma_v")%>% select(mean)
sigmav$mean
sigma_v
ypred <- fit_mcmc$summary("mu")
colMeans(matrix(ypred$mean,nrow = D) -Y)

