# Turn off scientific notation
options(scipen=999)

# Data
x <- c(9,1,2)
a <- c(0.000001, 0.000001, 0.000001)

# Posterior
posterior <- function(param) {
  lik <- sum(log(param^x))    # kernal of multinomial
  prior <- sum(log(param^a))  # kernal of Dirichlet
  return(lik + prior)
}


######## Metropolis algorithm ################

# proposal function
scale.factor=1000  # Scale variance of proposal function
proposalfunction <- function(param, scale.factor){
  # propose new values of param
  param <- rdirichlet(1,param*scale.factor)
  return(param)
}

# MCMC algorithm 
run_metropolis_MCMC <- function(startvalue, iterations, s = scale.factor){
  chain = array(dim = c(iterations+1,4))
  chain[1,1:3] = startvalue
  for (i in 1:iterations){
    # proposal function
    proposal = proposalfunction(chain[i,1:3], s)
    # ratio of likelihood density values
    plikelihood = posterior(proposal) - 
                  posterior(chain[i,1:3])
    ###### Ratio of proposal density values #######
    # Note: for some reason, this isn't needed,
    # even though a Dirichlet isn't perfectly symetrical
    # Verified by comparing to closed form solution
    #ppropos = log(ddirichlet(chain[i,1:3],proposal*s)) -
    #          log(ddirichlet(proposal,chain[i,1:3]*s))
    # probability of acceptence
    #probab = exp(plikelihood + ppropos)
    probab = exp(plikelihood)
    if(i %% (iterations/100) == 0) print(paste0((i/iterations)*100,"%"))
    if (runif(1) < probab){
      chain[i+1,1:3] = proposal
      chain[i+1,4] = posterior(proposal)
    }else{
      chain[i+1,1:3] = chain[i,1:3]
      chain[i+1,4] = posterior(chain[i,1:3])
    }
  }
  return(chain)
}

startvalue = x/sum(x)   # Initial values set to p of x
chain = run_metropolis_MCMC(startvalue, 100000)

burnIn = 1000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))
acceptance

### Plot Summary: #######################
# Plot of parameters
{
  par(mfrow = c(2,3))
  hist(chain[-(1:burnIn),1],nclass=30, , main="Posterior of p1", xlab="True value = red line" )
  abline(v = mean(chain[-(1:burnIn),1]), col="red")
  hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of p2", xlab="True value = red line")
  abline(v = mean(chain[-(1:burnIn),2]), col="red")
  hist(chain[-(1:burnIn),3],nclass=30, main="Posterior of p3", xlab="True value = red line")
  abline(v = mean(chain[-(1:burnIn),3]), col="red" )
  plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of p1", )
  plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of p2", )
  plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , main = "Chain values of p3", )
}

# Plot NPS
{
nps <- (chain[,1] - chain[,3]) * 100
par(mfrow = c(2,1))
hist(nps[-(1:burnIn)],nclass=30, , main="Posterior of p1", xlab="True value = red line" )
abline(v = mean(nps[-(1:burnIn)]), col="red")
plot(nps[-(1:burnIn)], type = "l" , main = "Chain values of NPS", )
mean(nps[-(1:burnIn)])
}


# Plot likelihood
{
  par(mfrow = c(2,1))
  hist(chain[-(1:burnIn),4],nclass=30, , main="Likelihood" )
  plot(chain[-(1:burnIn),4], type = "l" , main = "Chain values of likelihood", )
}
