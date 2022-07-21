#' ---
#' title: "Simple Regression Example"
#' author: Andrés Gutiérrez
#' date: "16th May 2021"
#' ---

#' # Setup
#+ message=FALSE

rm(list = ls())

library("rstan")
library("rstantools")
library(ggplot2)
library(tidyverse)
library(bayesplot)

rstan_options(auto_write = TRUE)
options(width = 90)


# data and simulated model ----------------------------------------------------------

#' # Generate data
N <- 30
x <- sort(rnorm(N))
a <- 2
b <- 3
sigma <- 1
y <- a + b * x + rnorm(N, 0, sigma)
sample_data <- list(N = N, x = x, y = y, sigma = sigma)

# STAN fit ----------------------------------------------------------------

#' # Draw from posterior distribution
# results='hide'
fit <- stan("5. modelos lineales/normal.stan", 
            data = sample_data, 
            verbose = TRUE)

#' ## Posterior summary and convergence diagnostics
print(fit, digits = 2)

# Plotting the MCMC -------------------------------------------------------
posterior <- as.array(fit)
dim(posterior)
dimnames(posterior)

thetapars <- c("a", "b", "sigma")
ypredpars <- dimnames(posterior)$parameters[4:33]

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

color_scheme_set("purple")
mcmc_dens(posterior, pars = thetapars)
mcmc_dens_overlay(posterior, pars = thetapars)

color_scheme_set("teal")
mcmc_violin(posterior, pars = thetapars, probs = c(0.1, 0.5, 0.9))

color_scheme_set("gray")
mcmc_scatter(posterior, pars = c("a", "b"), 
             size = 1.5, alpha = 0.5)

# requires hexbin package
if (requireNamespace("hexbin", quietly = TRUE)) {
  mcmc_hex(posterior, pars = c("a", "b"))
}


color_scheme_set("pink")
mcmc_pairs(posterior, pars = thetapars,
           off_diag_args = list(size = 1.5))

color_scheme_set("blue")
mcmc_trace(posterior, pars = thetapars)

color_scheme_set("mix-blue-red")
mcmc_trace(posterior, pars = thetapars, 
           facet_args = list(ncol = 1, strip.position = "left"))

mcmc_trace_highlight(posterior, pars = thetapars, highlight = 3)


# Visual MCMC diagnostics -------------------------------------------------

color_scheme_set("darkgray")
np <- nuts_params(fit)
mcmc_pairs(fit, pars = thetapars)

color_scheme_set("red")
mcmc_nuts_energy(np)

rhats <- rhat(fit)
print(rhats)

color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(rhats)

mcmc_acf(posterior, pars = thetapars, lags = 10)


# Graphical posterior predictive checks -----------------------------------

sims <- as.data.frame(fit)
yrep <- as.matrix(sims[, c(4: 33)])
rowsrandom <- sample(nrow(yrep), 200)
yrep2 <- yrep[rowsrandom, ]

color_scheme_set("brightblue")

ppc_dens_overlay(y, yrep2)
ppc_hist(y, yrep2)
ppc_ecdf_overlay(y, yrep2)

ppc_intervals(y, yrep2)
ppc_intervals(y, yrep)
ppc_ribbon(y, yrep2)
ppc_ribbon(y, yrep)
ppc_boxplot(y, yrep2)
ppc_scatter(y, yrep2)

prop_gzero <- function(x) mean(x >= 0)
prop_gzero(y) # check proportion of values greater tha zero in y
ppc_stat(y, yrep, stat = "prop_gzero")
ppc_stat(y, yrep, stat = "max")
ppc_stat(y, yrep, stat = "min")
ppc_stat(y, yrep, stat = "mean")
ppc_stat(y, yrep, stat = "sd")
ppc_stat_2d(y, yrep2)
ppc_stat(y, yrep2)


# Andres' own and inneficient graphs --------------------------------------

#' ## Plot the conditional mean function
sims <- as.data.frame(fit)
plot(x, y, las = 1, pch = 20)

for (s in 1:nrow(fit)) {
  with(sims, curve(a[s] + b[s] * x, lwd = 0.5, add = TRUE, col = "red"))
}

summary_fit <- get_posterior_mean(fit)
a.hat <- summary_fit[1,1]
b.hat <- summary_fit[2,1]
curve(a.hat + b.hat * x, add = TRUE, col = "green")
curve(a + b * x, add = TRUE, col = "blue")

sims1 <- sims %>% 
  select(-c(a, b, sigma, lp__)) 

colnames(sims1) <- as.character(1:N)

sims2 <- sims1  %>%
  pivot_longer(everything(),
               names_to = "iter",
               values_to = "y_pred") %>%
  mutate(k = sort(rep(1:nrow(sims), N)))

y.orig <- data.frame(
  y_pred = y,
  k = 0,
  i = 1:N
)

g1 <- sims2 %>% ggplot() +
  geom_density(
    aes(x = y_pred, 
        group = as.factor(k)),
    col = "grey"
  ) + 
  theme_classic() +
  theme(legend.position = "none")

g1 + geom_density(
  data = y.orig,
  aes(x = y_pred)
)

g2 <- sims2 %>% ggplot() +
  geom_boxplot(
    outlier.shape = NA,
    aes(y = y_pred, 
        x = iter)) + 
  theme_classic() +
  theme(legend.position = "none")

g2 + geom_point(
  data = y.orig,
  aes(x = reorder(i, i),
      y = y_pred),
  col = 2,
  size = 4
)

launch_shinystan(fit)
