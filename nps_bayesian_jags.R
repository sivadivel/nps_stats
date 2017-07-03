library(rjags)
library(coda)
# Housekeeping
{
  # Turn off scientific notation
  options(scipen=999)
  options(digits=4)
}

# Data
n <- c(9,0,2)
N <- sum(n)

# Model
model.text<-"
model {
  n[1:3] ~ dmulti(p[1:3],N)
  p[1:3] ~ ddirch(c(0.0001, 0.0001, 0.0001))
  nps <- 100 * (p[1] - p[3])
}
"
dat <- list('n'=n,'N'=N) # Data in list form for compiling model
its <- 10000 # Number of iterations

# MCMC
mcmc.sim <- function(mod, dat, its){
  # Create model
  model.sat.spec <- textConnection(mod)
  # Compile model
  sat.jags <- jags.model(model.sat.spec,
                       data=dat,
                       n.chains = 3,
                       n.adapt = 1000)
  # Sample from posterior
  samps.coda <- coda.samples(sat.jags,
                           c('p','nps'),
                           n.iter=its,
                           thin=1)
  return(samps.coda)
}

chains <- mcmc.sim(model.text, dat, its)

# Summary of results
results <- summary(chains)
# For some reason, only this column needs to be rounded
results[[2]][,5] <- round(results[[2]][,5], 4)
results

plot(chains[[1]][,1])

plot(chains[[1]][,c(2,4)])
