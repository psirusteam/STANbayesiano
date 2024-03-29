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
library(posterior)

rstan_options(auto_write = TRUE)
options(width = 90)

# data and simulated model ----------------------------------------------------------

#' # Generate data
n <- 500 # Tamaño de muestra
ns <- 100
nd <- round(runif(ns, 300,600)) # predicción
p <- 3  # Número de categorías
theta <- c(0.2, 0.45, 0.35) # parámetros de éxito 
sum(theta)
y <- t(rmultinom(n = n, size = 1, prob = theta))
#yagg <- colSums(y)
alpha = rep(0.5, 3)

sample_data <- list(n = n, p = p, y = y, alpha = alpha, ns = ns,nd = nd)

# STAN fit ----------------------------------------------------------------

#' # Draw from posterior distribution
#+ results='hide'
#' fit <- stan("2. multiparamétricos/Multinomial1.stan", 
#'             data = sample_data)
#' 
#' #' ## Posterior summary and convergence diagnostics
#' print(fit, digits = 2, pars = "theta")


fit <- cmdstan_model(stan_file = "2. multiparamétricos/Multinomial1.stan",
                     compile = TRUE)


fit_mcmc <- fit$sample(
  data = sample_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 0
)

fit_mcmc$summary()

draws_lin <- posterior::as_draws_df(fit_mcmc$draws())
draws_lin %>%
  as_draws_df() %>%
  as_tibble() %>%
  select(starts_with("ypred")) %>%
  apply(2, mean) %>%
  matrix(data = .,nrow = n,ncol = p) %>% 
  data.frame()  


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
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)

mcmc_areas(
  posterior, 
  pars = ypredpars,
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
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
mcmc_trace(posterior, pars = thetapars, 
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
yrep2 <- as.matrix(yrep[rowsrandom, ])

color_scheme_set("brightblue")

ppc_dens_overlay(y, yrep)
ppc_hist(y, yrep2)
ppc_ecdf_overlay(y, yrep)

prop_gzero <- function(x) mean(x == 0)
prop_gzero(y) # check proportion of values greater tha zero in y
prop_gones <- function(x) mean(x == 1)
prop_gones(y) # check proportion of values greater tha zero in y
ppc_stat(y, yrep, stat = "prop_gzero")
ppc_stat(y, yrep, stat = "prop_gones")
ppc_stat(y, yrep, stat = "mean")
ppc_stat(y, yrep, stat = "sd")

# Shiny checks ------------------------------------------------------------
library(shinystan)
launch_shinystan(fit)
