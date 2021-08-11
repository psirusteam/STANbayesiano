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

library(survival)
library(dplyr)
data(heart)

sobrevida <- heart %>%
  filter(transplant == 1) %>%
  mutate(tiempo = stop - start)

sample_data <- list(y = sobrevida$tiempo, 
                    n = nrow(sobrevida))

# STAN fit ----------------------------------------------------------------

fit <- stan("1. uniparamétricos/Exponencial.stan", 
            data = sample_data)

#' ## Posterior summary and convergence diagnostics
print(fit, digits = 4, pars = c("theta", "invtheta"))

bayesplot::mcmc_areas(fit, 
                      pars = c("theta", "invtheta"), 
                      prob = 0.95)

