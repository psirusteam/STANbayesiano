rm(list = ls())

library("rstan")
library("rstantools")
library(ggplot2)
library(tidyverse)
library(bayesplot)

rstan_options(auto_write = TRUE, 
              javascript = FALSE)
options(width = 90)
options(mc.cores = parallel::detectCores())

# data and simulated model ----------------------------------------------------------

Escuelas <- data.frame(
  row.names=c("A", "B", "C", "D", "E", "F", "G", "H"),
  efecto = c(28.39, 7.94, -2.75, 6.82,
             -0.64, 0.63, 18.01, 12.16))
Escuelas
hist(Escuelas$efecto)
mean(Escuelas$efecto)
sd(Escuelas$efecto)

sample_data <- list(y = Escuelas$efecto, 
                    n = nrow(Escuelas),
                    mu = mean(Escuelas$efecto))

# STAN fit ----------------------------------------------------------------

fit <- stan("1. uniparamétricos/NormalVarianza.stan", 
            data = sample_data)

#' ## Posterior summary and convergence diagnostics
print(fit, digits = 4, pars = c("tau", "tau2"))

bayesplot::mcmc_areas(fit, 
                      pars = c("tau"), 
                      prob = 0.95)

# Graphical posterior predictive checks -----------------------------------

sims <- as.data.frame(fit)
names(sims)
y <- sample_data$y
yrep <- as.matrix(sims[, 3:10])
rowsrandom <- sample(nrow(yrep), 20)
yrep2 <- as.matrix(yrep[rowsrandom, ])

color_scheme_set("brightblue")

ppc_dens_overlay(y, yrep)
ppc_dens_overlay(y, yrep2)

ppc_hist(y, yrep2)
ppc_ecdf_overlay(y, yrep)

prop_gzero <- function(x) mean(x >= 0)
prop_gzero(y) # check proportion of values greater tha zero in y
ppc_stat(y, yrep, stat = "prop_gzero")
ppc_stat(y, yrep, stat = "mean")
ppc_stat(y, yrep, stat = "sd")


