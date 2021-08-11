rm(list = ls())

library("rstan")
library("rstantools")
library(ggplot2)
library(tidyverse)
library(bayesplot)

rstan_options(auto_write = TRUE, 
              javascript = FALSE)
options(width = 90)

# data and simulated model ----------------------------------------------------------
n <- 12
y <- c(3.56, 3.36, 2.99, 2.71, 3.31, 3.68, 
       2.78, 2.95, 2.82, 3.45, 3.42, 3.15)

sample_data <- list(y = y, 
                    n = n)

# STAN fit ----------------------------------------------------------------

fit <- stan("1. uniparamétricos/NormalMedia.stan", 
            data = sample_data)

#' ## Posterior summary and convergence diagnostics
print(fit, digits = 4, pars = "theta")

bayesplot::mcmc_areas(fit, 
                      pars = c("theta", "invtheta"), 
                      prob = 0.95)

