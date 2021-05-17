library("bayesplot")
library("ggplot2")
library("rstanarm")

head(roaches) # see help("rstanarm-datasets")
roaches$roach100 <- roaches$roach1 / 100 # pre-treatment number of roaches (in 100s)

# using rstanarm's default priors. For details see the section on default
# weakly informative priors at https://mc-stan.org/rstanarm/articles/priors.html
fit_poisson <- stan_glm(
  y ~ roach100 + treatment + senior,
  offset = log(exposure2),
  family = poisson(link = "log"),
  data = roaches,
  seed = 1111,
  refresh = 0 # suppresses all output as of v2.18.1 of rstan
)

print(fit_poisson)

fit_nb <- update(fit_poisson, family = "neg_binomial_2")

print(fit_nb)

y <- roaches$y
y

yrep_poisson <- 
  posterior_predict(fit_poisson, draws = 500)
yrep_nb <- 
  posterior_predict(fit_nb, draws = 500)


dim(yrep_poisson)
dim(yrep_nb)

color_scheme_set("brightblue")
ppc_dens_overlay(y, yrep_poisson[1:50, ])
ppc_dens_overlay(y, yrep_poisson[1:50, ]) + xlim(0, 150)
ppc_hist(y, yrep_poisson[1:10, ])

ppc_dens_overlay(y, yrep_nb[1:50, ])+ xlim(0, 300)
ppc_hist(y, yrep_nb[1:5, ])
ppc_hist(y, yrep_nb[1:5, ], binwidth = 20) +
  coord_cartesian(xlim = c(-1, 300))

prop_zero <- function(x) mean(x == 0)
prop_zero(y) # check proportion of zeros in y

ppc_stat(y, yrep_poisson, stat = "prop_zero", binwidth = 0.005)
ppc_stat(y, yrep_nb, stat = "prop_zero")

ppc_stat(y, yrep_poisson, stat = "max")
ppc_stat(y, yrep_nb, stat = "max")
ppc_stat(y, yrep_nb, stat = "max", binwidth = 100) +
  coord_cartesian(xlim = c(-1, 5000))

available_ppc()
available_ppc(pattern = "_grouped")
ppc_stat_grouped(y, yrep_nb, group = roaches$treatment, stat = "prop_zero")
ppc_stat_grouped(y, yrep_poisson, group = roaches$treatment, stat = "prop_zero")




