library(ggplot2)
library(bayesplot)
library(rstanarm)

theme_set(bayesplot::theme_default())

data("womensrole", package = "HSAUR3")
womensrole$total <- 
  womensrole$agree + womensrole$disagree

womensrole_glm_1 <- 
  glm(cbind(agree, disagree) ~ education + gender,
      data = womensrole, family = binomial(link = "logit"))
round(coef(summary(womensrole_glm_1)), 3)

womensrole_bglm_1 <- 
  stan_glm(cbind(agree, disagree) ~ education + gender,
           data = womensrole,
           family = binomial(link = "logit")
           )

womensrole_bglm_1$algorithm

ci95 <- 
  posterior_interval(womensrole_bglm_1, 
                     prob = 0.95, 
                     pars = "education")
round(ci95, 2)

residuals(womensrole_bglm_1)

launch_shinystan(womensrole_bglm_1, ppd = FALSE)

rstan::get_stanmodel(womensrole_bglm_1$stanfit)    



flat_prior_test <- 
  stan_glm(mpg ~ wt, data = mtcars, prior = NULL)

rstan::get_stanmodel(flat_prior_test$stanfit)    


library(brms)
fit1 <- brm(count ~ zAge, 
            data = epilepsy)

summary(fit1) 
stancode(fit1)

make_stancode(rating ~ treat, 
              data = inhaler)



library("rstan") # observe startup messages
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

schools_dat <- list(J = 8, 
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

fit <- stan(file = 'schools.stan', data = schools_dat)

print(fit)
plot(fit)
pairs(fit, pars = c("mu", "tau", "lp__"))

la <- extract(fit, permuted = TRUE) # return a list of arrays 
mu <- la$mu 

### return an array of three dimensions: iterations, chains, parameters 
a <- extract(fit, permuted = FALSE) 

### use S3 functions on stanfit objects
a2 <- as.array(fit)
m <- as.matrix(fit)
d <- as.data.frame(fit)


