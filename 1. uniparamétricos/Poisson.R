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

#' # Generate data
y <- c(22, 9, 9, 20, 10, 14, 11, 14, 11, 11, 
       19, 12, 8, 9, 16, 8, 13, 8, 14, 12, 14, 
       11, 14, 13, 11, 14, 13, 11, 7, 12)

sample_data <- list(y = y, n = length(y))

# STAN fit ----------------------------------------------------------------

#' # Draw from posterior distribution
#+ results='hide'
fit1 <- stan("1. uniparamétricos/Poisson.stan", 
            data = sample_data)

#' ## Posterior summary and convergence diagnostics
print(fit1, digits = 5, pars = "theta")

bayesplot::mcmc_areas(fit1, pars = "theta", 
                      prob = 0.95)

