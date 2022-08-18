#' ---
#' title: "Simple Multinomial Example"
#' author: Andrés Gutiérrez
#' date: "22th May 2021"
#' ---

rm(list = ls())

library("rstan")
library("rstantools")
library(ggplot2)
library(tidyverse)
library(bayesplot)
library(cmdstanr)
library(brms)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
options(width = 90)

# data and simulated model ----------------------------------------------------------
# Definir parámetros
D <- 200 # número de dominios.
n <- round(runif(n = D, 200, 500)) # Tamaño de muestra por dominio
#n <- rep(1,D)
P <- 8 # número de indices.
## parametros de regresión.
K <- 2  # Número variables regresoras

beta <- t(matrix(rep(c(0.5,-0.2),P), ncol = P, nrow = K) )

X <- cbind(1, matrix(rnorm(D, 5, 3), 
                     D, (K - 1)))  #variables regresoras

### Descomposición de Cholesky para simulación de la
### correlación de los efectos aleatorios
C <- matrix(NA,P,P)
C <- 0.8^(abs(row(C)-col(C)))
L <- chol(C)
sigma <- diag(runif(P,min = 0.25, 0.8))
Lambda <- sigma %*% t(L)
Z <- rbind(rnorm(D), 
           rnorm(D),
           rnorm(D),
           rnorm(D),
           rnorm(D),
           rnorm(D),
           rnorm(D),
           rnorm(D))
u <- t(Lambda %*% Z)

cor(u)
sd(u[, 1])
sd(u[, 2])
mean(u[, 1])
mean(u[, 2])

### Calculando el denominador
lp1 <- crossprod(t(X), (beta[1, ])) + u[, 1]
lp2 <- crossprod(t(X), (beta[2, ])) + u[, 2]
lp3 <- crossprod(t(X), (beta[3, ])) + u[, 3]
lp4 <- crossprod(t(X), (beta[4, ])) + u[, 4]
lp5 <- crossprod(t(X), (beta[5, ])) + u[, 5]
lp6 <- crossprod(t(X), (beta[6, ])) + u[, 6]
lp7 <- crossprod(t(X), (beta[7, ])) + u[, 7]
lp8 <- crossprod(t(X), (beta[8, ])) + u[, 8]


logit1 <- exp(lp1)/(1+exp(lp1))
logit2 <- exp(lp2)/(1+exp(lp2))
logit3 <- exp(lp3)/(1+exp(lp3))
logit4 <- exp(lp4)/(1+exp(lp4))
logit5 <- exp(lp5)/(1+exp(lp5))
logit6 <- exp(lp6)/(1+exp(lp6))
logit7 <- exp(lp7)/(1+exp(lp7))
logit8 <- exp(lp8)/(1+exp(lp8))

y <- matrix(NA, D, K)
y <- as_tibble(y)

y <-  y %>%  mutate(
    V1 = map(1:D,function(ii){ rbinom(n = 1, size = n[ii], prob = logit1[ii])}),
    V2 = map(1:D,function(ii){ rbinom(n = 1, size = n[ii], prob = logit2[ii])}),
    V3 = map(1:D,function(ii){ rbinom(n = 1, size = n[ii], prob = logit3[ii])}),
    V4 = map(1:D,function(ii){ rbinom(n = 1, size = n[ii], prob = logit4[ii])}),
    V5 = map(1:D,function(ii){ rbinom(n = 1, size = n[ii], prob = logit5[ii])}),
    V6 = map(1:D,function(ii){ rbinom(n = 1, size = n[ii], prob = logit6[ii])}),
    V7 = map(1:D,function(ii){ rbinom(n = 1, size = n[ii], prob = logit7[ii])}),
    V8 = map(1:D,function(ii){ rbinom(n = 1, size = n[ii], prob = logit8[ii])})) %>% 
unnest(cols = c(V1, V2, V3, V4, V5, V6, V7, V8))
y$nd <- n

sample_data <- list(D = D,
                    P = P,
                    K = K,
                    Kz = D,
                    y = y[,-9] %>% as.matrix(),
                    X = X,
                    Z = diag(1,D),
                    nd = n)
# STAN fit ----------------------------------------------------------------

#' # Draw from posterior distribution
#+ results='hide'

# 
#  fit <- stan("Modelos/Multinomial7.stan",            
#              data = sample_data)

fit <-
  cmdstan_model(stan_file = "2. multiparamétricos/MultiBinomialCorrelacionado.stan",
                compile = TRUE)


fit_mcmc <- fit$sample(
  data = sample_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4
)
 fit_mcmc$save_object(
   file = "fit1.rds")


fit_mcmc$print("sigma_u")
# fit_mcmc$print("sigma_u2")
fit_mcmc$print("beta")
fit_mcmc$print("theta")
#' ## Posterior summary and convergence diagnostics
#' 

draws <- as_draws_df(fit_mcmc$draws())

betas <- draws %>% 
  select(matches("beta")) %>% 
  colMeans() %>% 
  matrix(., ncol = K,
         nrow = P,
         byrow = FALSE)

omega <- draws %>% 
  select(matches("^Omega")) %>% 
  colMeans() %>% 
  matrix(., ncol = P,
         nrow = P,
         byrow = FALSE)
cor(u)

draws %>% 
  select(matches("^sdcomprobar")) %>% 
  colMeans() %>% 
  matrix(., ncol = P,
         nrow = 1,
         byrow = FALSE)
apply(u, 2, sd)


efecto_aleatorio <- draws %>% 
  select(matches("^u")) %>% 
  colMeans() %>% 
  matrix(., ncol = P,
         nrow = D,
         byrow = FALSE)
apply(u, 2, sd)
apply(efecto_aleatorio, 2, sd)

thetas <- draws %>% 
  select(matches("^theta")) %>% 
  colMeans() %>% 
  matrix(., ncol = P,
         nrow = D,
         byrow = FALSE)
par(mfrow = c(2,4))
plot(thetas[,1],logit1)
abline(a = 0, b =1, col = "red")
plot(thetas[,2],logit2)
abline(a = 0, b =1, col = "red")
plot(thetas[,3],logit3)
abline(a = 0, b =1, col = "red")
plot(thetas[,4],logit4)
abline(a = 0, b =1, col = "red")
plot(thetas[,5],logit5)
abline(a = 0, b =1, col = "red")
plot(thetas[,6],logit6)
abline(a = 0, b =1, col = "red")
plot(thetas[,7],logit7)
abline(a = 0, b =1, col = "red")
plot(thetas[,8],logit8)
abline(a = 0, b =1, col = "red")


print(fit, digits = 2, pars = "theta")
print(fit, digits = 2, pars = "beta")
print(fit, digits = 2, pars = "sigma_u")
print(fit, digits = 2, pars = "Omega")
print(fit, digits = 2, pars = "u")
print(fit, digits = 2, pars = "z_u")
print(fit, digits = 2, pars = "tau_u")


# Plotting the MCMC -------------------------------------------------------
posterior <- as.array(fit)
dim(posterior)
dimnames(posterior)

pos <- grep(x = dimnames(posterior)$parameters,
      pattern = "beta")

betapars <- dimnames(posterior)$parameters[pos]

color_scheme_set("red")
plot(fit, pars = betapars)
plot(fit, pars = ypredpars)

mcmc_intervals(posterior, pars = betapars)

mcmc_areas(
  posterior,
  pars = betapars,
  prob = 0.8,
  # 80% intervals
  prob_outer = 0.99,
  # 99%
  point_est = "mean"
)

color_scheme_set("green")
mcmc_hist(posterior, pars = betapars)

color_scheme_set("brightblue")
mcmc_hist_by_chain(posterior, pars = betapars)


color_scheme_set("purple")
mcmc_dens(posterior, pars = betapars)
mcmc_dens_overlay(posterior, pars = betapars)

color_scheme_set("teal")
mcmc_violin(posterior,
            pars = betapars, probs = c(0.1, 0.5, 0.9))

color_scheme_set("blue")
mcmc_trace(posterior, pars = betapars)
mcmc_trace(posterior, pars = ypredpars)

color_scheme_set("mix-blue-red")
mcmc_trace(posterior,
           pars = betapars,
           facet_args = list(ncol = 1, strip.position = "left"))

mcmc_trace_highlight(posterior, pars = betapars,
                     highlight = 3)

# Visual MCMC diagnostics -------------------------------------------------

color_scheme_set("red")

rhats <- rhat(fit)
print(rhats)

color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(rhats)

mcmc_acf(posterior, pars = betapars, lags = 10)

# Graphical posterior predictive checks -----------------------------------

sims <- as.data.frame(fit)
yrep <- as.matrix(sims[, 2:31])
rowsrandom <- sample(nrow(yrep), 20)
yrep2 <- as.matrix(yrep[rowsrandom,])

color_scheme_set("brightblue")

ppc_dens_overlay(y, yrep)
ppc_hist(y, yrep2)
ppc_ecdf_overlay(y, yrep)

prop_gzero <- function(x)
  mean(x == 0)
prop_gzero(y) # check proportion of values greater tha zero in y
prop_gones <- function(x)
  mean(x == 1)
prop_gones(y) # check proportion of values greater tha zero in y
ppc_stat(y, yrep, stat = "prop_gzero")
ppc_stat(y, yrep, stat = "prop_gones")
ppc_stat(y, yrep, stat = "mean")
ppc_stat(y, yrep, stat = "sd")

# Shiny checks ------------------------------------------------------------
library(shinystan)
launch_shinystan(fit)
