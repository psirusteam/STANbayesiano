# BUGS program for Binomial-logistic regression model with AR(1) dependence.
# Data is for m areas and T time points. Regression model has distinct 
# intercepts and other regrerssion coefficients for each year.

# Note: Data (neff, yeff, IRSpr, IRSnf, FS), for all m areas and 
# all T years in m x T arrays.

#---------------------------------------------------------------------
# Model specification section

  model {

    for (i in 1:m) {

# First year of data

        y[i,1] ~ dbin(p[i,1], n[i,1])

        w[i,1] ~ dnorm(mu[i,1],tau1)
        logit(p[i,1]) <- w[i,1]

        mu[i,1] <- b0[1] + beta1[1]*(IRSpr[i,1] - mean(IRSpr[,1]))
                         + beta2[1]*(IRSnf[i,1] - mean(IRSnf[,1]))
                         + beta3[1]*(FS[i,1] - mean(FS[,1]))

# Subsequent years of data

      for (j in 2:T) {

        y[i,j] ~ dbin(p[i,j], n[i,j])

        w[i,j] ~ dnorm(condmu[i,j],tau)
        logit(p[i,j]) <- w[i,j]

        mu[i,j] <- b0[j] + beta1[j]*(IRSpr[i,j] - mean(IRSpr[,j]))
                         + beta2[j]*(IRSnf[i,j] - mean(IRSnf[,j]))
                         + beta3[j]*(FS[i,j] - mean(FS[,j]))

        condmu[i,j] <- mu[i,j] + phi*(w[i,j-1] - mu[i,j-1])
      }
    }


# Set tau1, precision for first time point.

    tau1 <- tau*(1-phi*phi)


# Prior specifications for regression coefficients

    b0[1:T] ~ dmnorm(mub0[],Rb0[,])
    beta1[1:T] ~ dmnorm(mubeta[],Rbeta[,])
    beta2[1:T] ~ dmnorm(mubeta[],Rbeta[,])
    beta3[1:T] ~ dmnorm(mubeta[],Rbeta[,])

# Uniform prior for the AR(1) parameter phi

    phi ~ dunif(-1,1)
#   phi ~ dunif(0,1) 

# Alternative prior specifications for the model error variance. 

# Inverse Gamma prior (Gamma prior on tau)
#    tau ~ dgamma(.01,.01)
#    sigma <- 1/sqrt(tau)
#    s2 <- 1/tau

# Uniform prior on model error variance s2.
#    s2 ~ dunif(0,1)
#    tau <- 1 / s2
#    sigma <- sqrt(s2)

# Uniform prior on model error std. deviation sigma.
   sigma ~ dunif(0,3)
   s2 <- sigma*sigma
   tau <- 1 / s2

  }

#---------------------------------------------------------------------
# Data input section

# Set number of areas and number of time points.

#  list(m = 51) 
  list(m = 51,T = 6)   

# Set up zero vector and identity matrix for prior specification of b0
  list(mub0 = c(0,0,0,0,0,0),
  Rb0 = structure( .Data = c(.001,  0,   0,   0,   0,   0,
                               0, .001,  0,   0,   0,   0,
                               0,   0, .001,  0,   0,   0,
                               0,   0,   0, .001,  0,   0,
                               0,   0,   0,   0, .001,  0,
                               0,   0,   0,   0,   0, .001),
                           .Dim = c(6,6)))

# Set up zero vector and identity matrix for prior specification of beta
  list(mubeta = c(0,0,0,0,0,0),
  Rbeta = structure( .Data = c(.001,  0,   0,   0,   0,   0,
                                 0, .001,  0,   0,   0,   0,
                                 0,   0, .001,  0,   0,   0,
                                 0,   0,   0, .001,  0,   0,
                                 0,   0,   0,   0, .001,  0,
                                 0,   0,   0,   0,   0, .001),
                           .Dim = c(6,6)))

#---------------------------------------------------------------------
# Compile model (Model -> Specification Tool -> compile)


#---------------------------------------------------------------------
# Initial Values Section

# Set initial values for fixed effect parameters. Load these once for
# each MCMC chain to be run. (Model -> Specification Tool -> load inits).
# Use gen inits to generate initial values for random effects parameters.

#  list(b0 = c(0,0,0,0,0,0), beta1 = c(0,0,0,0,0,0), beta2 = c(0,0,0,0,0,0), 
#     beta3 = c(0,0,0,0,0,0), tau = 1, phi=.6)
#  list(b0 = c(-1.2,-1.2,-1.2,-1.2,-1.2,-1.2), beta1 = c(.9,.9,.9,.9.,9.,9), 
#     beta2 = c(.16,.16,.16,.16,.16,.16), beta3=c(.36,.36,.36,.36,.36,.36), tau = 30, phi=0)
#  list(b0 = c(-1.5,-1.5,-1.5,-1.5,-1.5,-1.5), beta1 = c(1.5,1.5,1.5,1.5,1.5,1.5), 
#     beta2 = c(.5,.5,.5,.5,.5,.5), beta3 = c(1,1,1,1,1,1), tau = 100, phi=.3)

#  list(b0 = c(0,0,0,0,0,0), beta1 = c(0,0,0,0,0,0), 
#     beta2 = c(0,0,0,0,0,0), beta3 = c(0,0,0,0,0,0), s2 = 1,phi=.6)
#  list(b0 = c(-1.2,-1.2,-1.2,-1.2,-1.2,-1.2), beta1 = c(.9,.9,.9,.9.,9.,9), 
#     beta2 = c(.16,.16,.16,.16,.16,.16), beta3=c(.36,.36,.36,.36,.36,.36), s2 = .0324, phi=0)
#  list(b0 = c(-1.5,-1.5,-1.5,-1.5,-1.5,-1.5), beta1 = c(1.5,1.5,1.5,1.5,1.5,1.5), 
#     beta2 = c(.5,.5,.5,.5,.5,.5), beta3 = c(1,1,1,1,1,1), s2 = .01, phi=.3)

  list(b0 = c(0,0,0,0,0,0), beta1 = c(0,0,0,0,0,0), 
     beta2 = c(0,0,0,0,0,0), beta3 = c(0,0,0,0,0,0), sigma = 1, phi=.6)
  list(b0 = c(-1.2,-1.2,-1.2,-1.2,-1.2,-1.2), beta1 = c(.9,.9,.9,.9.,9.,9), 
     beta2 = c(.16,.16,.16,.16,.16,.16), beta3=c(.36,.36,.36,.36,.36,.36), sigma = .18, phi=0)
  list(b0 = c(-1.5,-1.5,-1.5,-1.5,-1.5,-1.5), beta1 = c(1.5,1.5,1.5,1.5,1.5,1.5), 
     beta2 = c(.5,.5,.5,.5,.5,.5), beta3 = c(1,1,1,1,1,1), sigma = .05, phi=.3)
