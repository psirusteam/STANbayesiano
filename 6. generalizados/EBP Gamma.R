rm(list = ls())

library("rstan")
library("rstantools")
library(ggplot2)
library(tidyverse)
library(bayesplot)
library(magrittr)

rstan_options(auto_write = TRUE)
options(width = 90)
set.seed(2018)

N <- 1000

depto <- c(1:23)
depto <- sample(depto,
                N,
                replace = TRUE,
                prob = c(rep(0.1, 3), rep(0.05, 10), rep(0.02, 10)))
prop.table(table(depto))

area <- c(1:2)
area <- sample(area,
               N,
               replace = TRUE,
               prob = c(0.3, 0.7))
prop.table(table(area))

etnia <- c(1:3)
etnia <- sample(etnia,
                N,
                replace = TRUE,
                prob = c(0.1, 0.2, 0.70))
prop.table(table(etnia))

edad <- c(1:5)
edad <-
  sample(edad,
         N,
         replace = TRUE,
         prob = c(0.22, 0.25, 0.20, 0.22, 0.11))
prop.table(table(edad))

sexo <- c(1:2)
sexo <- sample(sexo,
               N,
               replace = TRUE,
               prob = c(0.51, 0.49))
prop.table(table(sexo))

anoest <- c(1:4)
anoest <-
  sample(anoest,
         N,
         replace = TRUE,
         prob = c(0.09, 0.28, 0.4, 0.23))
prop.table(table(anoest))

beta_0 <- 0.8
beta_depto <- data.frame(depto = as.factor(c(1:23)),
                         beta_depto = seq(-1.1, 0, by = 0.05))
beta_area <- data.frame(area = as.factor(c(1:2)),
                        beta_area = c(-0.050, 0))
beta_etnia <- data.frame(etnia = as.factor(c(1:3)),
                         beta_etnia = seq(-0.1, 0, by = 0.05))
beta_edad <- data.frame(edad = as.factor(c(1:5)),
                        beta_edad = seq(-0.2, 0, by = 0.05))
beta_sexo <- data.frame(sexo = as.factor(c(1:2)),
                        beta_sexo = seq(-0.050, 0, by = 0.05))
beta_anoest <- data.frame(anoest = as.factor(c(1:4)),
                          beta_anoest = c(-0.100, -0.075, -0.025, 0))

Xfactor <- data.frame(
  depto = as.factor(depto),
  area = as.factor(area),
  etnia = as.factor(etnia),
  edad = as.factor(edad),
  sexo = as.factor(sexo),
  anoest = as.factor(anoest)
)

Xdummy <-
  fastDummies::dummy_cols(Xfactor)

DGP <- left_join(Xdummy, beta_depto, by = "depto") %>%
  left_join(beta_area, by = "area") %>%
  left_join(beta_etnia, by = "etnia") %>%
  left_join(beta_edad, by = "edad") %>%
  left_join(beta_sexo, by = "sexo") %>%
  left_join(beta_anoest, by = "anoest") %>%
  mutate(
    LP = beta_0 + beta_depto + beta_area + beta_etnia + 
      beta_edad + beta_sexo + beta_anoest,
    theta = exp(LP), # función de vínculo es la identidad
    phi = 0.5,
    a = 1 / phi,
    b = a / (theta),
    y = rgamma(n = N, shape = a, rate = b)
  ) %>% 
  select(-depto, -area, -etnia, -edad, -sexo, -anoest)

hist(DGP$LP)
summary(DGP$LP)

hist(DGP$theta)
summary(DGP$theta)

hist(DGP$y)
summary(DGP$y)
sd(DGP$y)

Censofactor <- data.frame(ingreso = DGP$y, Xfactor)
Xdummy %<>% select(-depto, -area, -etnia, -edad, -sexo, -anoest,
                   -depto_23, -area_2, -etnia_3, 
                   -edad_5, -sexo_2, -anoest_4)
Censodummy <- data.frame(ingreso = DGP$y, b0 = 1, Xdummy)
# Censodummy <- data.frame(LP = DGP$LP, Xdummy)
# lm(LP ~ 1 + ., data =  Censodummy)


Censofactor %>%
  group_by(depto) %>%
  summarise(ingreso.medio = mean(ingreso)) %>%
  View()

Censofactor %>%
  group_by(area) %>%
  summarise(ingreso.medio = mean(ingreso)) %>%
  View()

Censofactor %>%
  group_by(etnia) %>%
  summarise(ingreso.medio = mean(ingreso)) %>%
  View()

Censofactor %>%
  group_by(edad) %>%
  summarise(ingreso.medio = mean(ingreso)) %>%
  View()

Censofactor %>%
  group_by(sexo) %>%
  summarise(ingreso.medio = mean(ingreso)) %>%
  View()

Censofactor %>%
  group_by(anoest) %>%
  summarise(ingreso.medio = mean(ingreso)) %>%
  View()

saveRDS(Censofactor, "2. Rcodes/STAN/Datafactor.rds")
saveRDS(Censodummy, "2. Rcodes/STAN/Datadummy.rds")

n <- 500
N <- N

encuestadummy <- Censodummy %>%
  sample_n(n)

Xsam <- encuestadummy %>% select(-ingreso)
Xpred <- Censodummy %>% select(-ingreso)

gamma_data <- list(
  n = n,
  # tamaño de muestra
  N = N,
  # tamaño poblacional
  p = ncol(Xsam),
  # número de predictores
  Xsam = Xsam,
  # matriz de diseño para el modelo de estimación
  Xpred = Xpred,
  # matriz para el modelo de predicción
  y = encuestadummy$ingreso # variable respuesta
)

# STAN fit ----------------------------------------------------------------

#' # Draw from posterior distribution
#+ results='hide'
fit <- stan(
  "2. Rcodes/STAN/EBP_gamma.stan",
  data = gamma_data,
  verbose = TRUE,
  cores = 8
)

#' ## Posterior summary and convergence diagnostics
print(fit, digits = 4, pars = "phi") 
unique(DGP$phi)
betamean = as.array(summary(fit, pars = "beta"))$summary[, 1]
betamean
sum(betamean)

print(fit, digits = 4, pars = "beta") 
print(fit, digits = 4, pars = "ysam") 
print(fit, digits = 4, pars = "ypred") 

# SAM - Graphical posterior predictive checks -----------------------------------
posterior <- as.data.frame(fit)
y <- gamma_data$y

id <- str_detect(names(posterior), "ysam")
ysam <- posterior[, id]
rowsrandom <- sample(nrow(ysam), 100)
ypredsam <- as.matrix(ysam[rowsrandom, ])

color_scheme_set("brightblue")
ppc_dens_overlay(y = y, ypredsam)

# PRED - Graphical posterior predictive checks -----------------------------------
posterior <- as.data.frame(fit)
y <- Censofactor$ingreso

id <- str_detect(names(posterior), "ypred")
ypred <- posterior[, id]
rowsrandom <- sample(nrow(ypred), 100)
ypredcenso <- as.matrix(ypred[rowsrandom, ])

color_scheme_set("brightblue")
ppc_dens_overlay(y = y, ypredcenso)
