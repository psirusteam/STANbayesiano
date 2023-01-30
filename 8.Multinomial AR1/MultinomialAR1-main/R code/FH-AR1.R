# BUGS program for Fay-Herriot linear model with AR(1) dependence.
# Data is for m areas and T time points. Regression model has distinct 
# intercepts for each year, but the other regression coefficients are 
# common across years.

# Note: Program is set to use m = 51 areas and T = 6 years of data. 
# These can be changed in the data input section. Changing either m 
# or T will require changing dimensions of other arrays specified 
# in the data input section.

# Note: Program is set up for using 6 years of state CPS 5-17 total 
# poverty rate data. Model does not use previous census estimates or 
# residuals.

#---------------------------------------------------------------------
# Model specification section

  model {

   for (i in 1:m) {

# First year of data

     y[i,1] ~ dnorm(w[i,1],tauy[i,1])

       w[i,1] ~ dnorm(mu[i,1],tau1)

       mu[i,1] <- b0[1] + beta[1]*(IRSpr[i,1] - mean(IRSpr[,1]))
                        + beta[2]*(IRSnf[i,1] - mean(IRSnf[,1]))
                        + beta[3]*(FS[i,1] - mean(FS[,1]))

       tauy[i,1] <- 1/v[i,1]

# Subsequent years of data

     for (j in 2:T) {

       y[i,j] ~ dnorm(w[i,j],tauy[i,j])

         w[i,j] ~ dnorm(condmu[i,j],tau)

         mu[i,j] <- b0[j] + beta[1]*(IRSpr[i,j] - mean(IRSpr[,j]))
                          + beta[2]*(IRSnf[i,j] - mean(IRSnf[,j]))
                          + beta[3]*(FS[i,j] - mean(FS[,j]))

         condmu[i,j] <- mu[i,j] + phi*(w[i,j-1] - mu[i,j-1])

         tauy[i,j] <- 1/v[i,j]
      }
    }


# Set tau1, the unconditional precision of w[i,j], which is determined 
# from tau and phi by the stationarity assumption. It is needed in the 
# model for w[i,1] (first year). Here tau, the innovation precision is 
# the one precision parameter to be estimated. Alternatively, we can let 
# tau1 be the one precision parameter to be estimated and let tau be 
# determined from tau1 and phi.

    tau1 <- tau*(1-phi*phi)
#    tau <- tau1/(1-phi*phi)


# Prior specifications for intercept terms b0[1:T] and beta[1:3]

    b0[1:T] ~ dmnorm(mub0[],Rb0[,])
    beta[1:3] ~ dmnorm(mubeta[],Rbeta[,])

# Uniform prior for the AR(1) parameter phi

    phi ~ dunif(-1,1)
#   phi ~ dunif(0,1) 

# Alternative prior specifications for the model error variance. 

# Inverse Gamma prior (Gamma prior on tau or tau1)

#   tau ~ dgamma(.01,.01)
#   sigma <- 1/sqrt(tau)
#   s2 <- 1/tau

#   tau1 ~ dgamma(.01,.01)
#   sigma <- 1/sqrt(tau1)
#   s2 <- 1/tau1

# Uniform prior on model error variance s2.
   s2 ~ dunif(0,25)
   tau <- 1 / s2
   sigma <- sqrt(s2)

# Uniform prior on model error std. deviation sigma.
#   sigma ~ dunif(0,5)
#   s2 <- sigma*sigma
#   tau <- 1 / s2

  }

#---------------------------------------------------------------------
# Data input section

# Set number of areas and number of time points.

#  list(m = 51) 
  list(m = 51,T = 6)   

# Set up zero vector and identity matrix for prior specification of b0. 
# Make sure the dimensions here are T x 1 and T x T.

  list(mub0 = c(0,0,0,0,0,0),
  Rb0 = structure( .Data = c(.001,  0,   0,   0,   0,   0,
                               0, .001,  0,   0,   0,   0,
                               0,   0, .001,  0,   0,   0,
                               0,   0,   0, .001,  0,   0,
                               0,   0,   0,   0, .001,  0,
                               0,   0,   0,   0,   0, .001),
                           .Dim = c(6,6)))

# Set up zero vector and identity matrix for prior specification of beta
  list(mubeta = c(0,0,0),
  Rbeta = structure( .Data = c(.001,  0,   0,
                                 0, .001,  0,
                                 0,   0, .001),
                           .Dim = c(3,3)))


#---------------------------------------------------------------------
# Compile model (Model -> Specification Tool -> compile)


#---------------------------------------------------------------------
# Initial Values Section

# Set initial values for fixed effect parameters. Load these once for
# each MCMC chain to be run. (Model -> Specification Tool -> load inits).
# Use gen inits to generate initial values for random effects parameters.

#  list(b0=c(20,20,20,20,20,20), beta = c(0,0,0), tau = 1, phi=.5)
#  list(b0=c(15,15,15,15,15,15), beta = c(.9,.16,.36), tau = .5, phi=0)
#  list(b0=c(25,25,25,25,25,25), beta = c(1.5,.5,1), tau = 2.0, phi=.3)

#  list(b0=c(20,20,20,20,20,20), beta = c(0,0,0), tau1 = .75, phi=.5)
#  list(b0=c(15,15,15,15,15,15), beta = c(.9,.16,.36), tau1 = .5, phi=0)
#  list(b0=c(25,25,25,25,25,25), beta = c(1.5,.5,1), tau1 = 1.82, phi=.3)

  list(b0=c(20,20,20,20,20,20), beta = c(0,0,0), s2 = 1,phi=.5)
  list(b0=c(15,15,15,15,15,15), beta = c(.9,.16,.36), s2 = 2.0, phi=0)
  list(b0=c(25,25,25,25,25,25), beta = c(1.5,.5,1), s2 = .5, phi=.3)

#  list(b0=c(20,20,20,20,20,20), beta = c(0,0,0), sigma = 1, phi=.5)
#  list(b0=c(15,15,15,15,15,15), beta = c(.9,.16,.36), sigma = 1.41, phi=0)
#  list(b0=c(25,25,25,25,25,25), beta = c(1.5,.5,1), sigma = .71, phi=.3)

