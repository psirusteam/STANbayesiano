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
y <- 1106
k <- 5
sample_data <- list(k = k, y = y)

# STAN fit ----------------------------------------------------------------

#' # Draw from posterior distribution
#+ results='hide'
fit <- stan("1. uniparamétricos/BinomialNegativa.stan", 
            data = sample_data)

#' ## Posterior summary and convergence diagnostics
print(fit, digits = 5, pars = "theta")


# Muestra de datos --------------------------------------------------------

y <- c(1001, 978, 999, 860, 1155, 585, 1030, 
       960, 1002, 763, 1036, 779, 1158, 1017, 
       888, 977, 1256, 1349, 1047, 1088, 649, 
       765, 699, 1042, 1212, 660, 671, 835, 
       997, 1146, 1016)
k <- c(4, 6, 5, 4, 4, 6, 3, 5, 6, 7, 5, 5, 4, 
       5, 6, 4, 6, 6, 5, 5, 3, 4, 5, 4, 5, 5, 
       5, 6, 5, 4, 5)
sample_data <- list(k = k, y = y, n = length(y))

fit <- stan("1. uniparamétricos/BinomialNegativa2.stan", 
            data = sample_data)

#' ## Posterior summary and convergence diagnostics
print(fit, digits = 5, pars = "theta")
