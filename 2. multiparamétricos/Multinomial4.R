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

rstan_options(auto_write = TRUE)
options(width = 90)

# data and simulated model ----------------------------------------------------------
# Definir parámetros
D <- 10 # número de dominios.
n <- round(runif(n = D, 200, 500)) # Tamaño de muestra por dominio
P <- 3  # Número de categorías

## parametros de regresión.
K <- 2 # número de variables regresoras.
beta <- matrix(c(-1, -2, 1, 1), P-1, K) # coefientes del modelo.
X <- cbind(1, matrix(rnorm(D, 5, 3), 
                     D, (K - 1)))  #variables regresoras
### Calculando el denominador
Denominador <- 1 + exp((crossprod(t(X), (beta[1, ])))) +
  exp((crossprod(t(X), (beta[2, ]))))

# Calculando la matrix y
theta <- cbind(1 / Denominador,
               exp((crossprod(t(X), (beta[1, ])))) / Denominador,
               exp((crossprod(t(X), (beta[2, ])))) / Denominador)

theta
rowSums(theta)
y <- matrix(NA, D, P)
for (ii in 1:D) {
  y[ii, ] <- t(rmultinom(n = 1, size = n[ii], prob = theta[ii, ]))
}
rowSums(y)

sample_data <- list(D = D,
                    P = P,
                    K = K,
                    y = y,
                    X = X)


# STAN fit ----------------------------------------------------------------

#' # Draw from posterior distribution
#+ results='hide'
fit <- stan("2. multiparamétricos/Multinomial4.stan",
            data = sample_data)

# fit <-
#   cmdstan_model(stan_file = "2. multiparamétricos/Multinomial4.stan",
#                 compile = TRUE)
# 
# 
# fit_mcmc <- fit$sample(
#   data = sample_data,
#   seed = 123,
#   chains = 4,
#   parallel_chains = 4
# )
#fit_mcmc$print("theta")

#' ## Posterior summary and convergence diagnostics
print(fit, digits = 2, pars = "theta")
print(fit, digits = 2, pars = "beta")

# Plotting the MCMC -------------------------------------------------------
posterior <- as.array(fit)
dim(posterior)
dimnames(posterior)

thetapars <- c("theta")
ypredpars <- dimnames(posterior)$parameters[2:31]

color_scheme_set("red")
plot(fit, pars = thetapars)
plot(fit, pars = ypredpars)

mcmc_intervals(posterior, pars = thetapars)
mcmc_intervals(posterior, pars = ypredpars)

mcmc_areas(
  posterior,
  pars = thetapars,
  prob = 0.8,
  # 80% intervals
  prob_outer = 0.99,
  # 99%
  point_est = "mean"
)

mcmc_areas(
  posterior,
  pars = ypredpars,
  prob = 0.8,
  # 80% intervals
  prob_outer = 0.99,
  # 99%
  point_est = "mean"
)

color_scheme_set("green")
mcmc_hist(posterior, pars = thetapars)
mcmc_hist(posterior, pars = ypredpars)

color_scheme_set("brightblue")
mcmc_hist_by_chain(posterior, pars = thetapars)
mcmc_hist_by_chain(posterior, pars = ypredpars)

color_scheme_set("purple")
mcmc_dens(posterior, pars = thetapars)
mcmc_dens_overlay(posterior, pars = thetapars)

color_scheme_set("teal")
mcmc_violin(posterior,
            pars = thetapars, probs = c(0.1, 0.5, 0.9))

color_scheme_set("blue")
mcmc_trace(posterior, pars = thetapars)
mcmc_trace(posterior, pars = ypredpars)

color_scheme_set("mix-blue-red")
mcmc_trace(posterior,
           pars = thetapars,
           facet_args = list(ncol = 1, strip.position = "left"))

mcmc_trace_highlight(posterior, pars = thetapars,
                     highlight = 3)

# Visual MCMC diagnostics -------------------------------------------------

color_scheme_set("red")

rhats <- rhat(fit)
print(rhats)

color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(rhats)

mcmc_acf(posterior, pars = thetapars, lags = 10)

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
