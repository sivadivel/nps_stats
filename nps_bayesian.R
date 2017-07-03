library(MCMCpack)

##### CREDIBLE INTERVALS ON ONE NPS SCORE #################

# Assuming a multinomial distribution for the number of
# detractors, neutrals and promoters, and a Dirichlet prior for the parameters
# of the multinomial, the posterior distribution of the counts of people in 
# each category is distributed Dirichlet with parameters =  alpha(i) + x(i)

# Take 100k samples from posterior and convert to NPS score
# np = # of promoters; nq = # of neutrals; nd = # of detractors
nps.sim <- function(np,nq,nd){
  # Parameters of minimally informative Dirichlet prior
  alpha <- c(0.000001, 0.000001, 0.000001)
  
  # Parameters of Dirichlet posteriors 
  theta <- c(np, nq, nd) + alpha
  
  ns <- 100000 # Number of simulations 
  
  # n random draws from Dirichlet(theta)
  sim <- rdirichlet(ns, theta)
  
  # Convert to NPS scores
  nps <- 100 * (sim[,1] - sim[,3])
  return(nps)
}

# Calculate point estimate and credible intervals
nps.results <- function(sim){
  # Point estimate and credible intervals
  point <- mean(sim)
  lb <- quantile(sim, 0.025)
  ub <- quantile(sim, 0.975)
  
  out <- list(lb, point, ub)
  names(out) <- c('Lower Bound', 'Point Estimate', 'Upper Bound')
  return(out)
}

nps <- nps.sim(9,0,2)
results <-nps.results(nps)
results

# Plot posterior
hist(nps, breaks=100)
{abline(v=results[1], col='red', lwd=3)
abline(v=mean(nps), col='red', lwd=3)
abline(v=results[3], col='red', lwd=3)}


##### CREDIBLE INTERVALS ON THE DIFFERENCE BETWEEN TWO NPS SCORES #######

# Take 100k samples from posterior and convert to NPS score
# and convert to difference in NPS scores
nps.sim.2test <- function(np1,nq1,nd1,np2,nq2,nd2){
  # Parameters of minimally informative Dirichlet prior
  alpha1 <- c(0.00001, 0.00001, 0.00001)
  alpha2 <- c(0.00001, 0.00001, 0.00001)
  
  # Parameters of Dirichlet posteriors 
  theta1 <- c(np1, nq1, nd1) + alpha1
  theta2 <- c(np2, nq2, nd2) + alpha2
  
  ns <- 100000 # Number of simulations 
  
  # n random draws from Dirichlet(theta)
  sim1 <- rdirichlet(ns, theta1)
  sim2 <- rdirichlet(ns, theta2)
  
  # Convert to NPS scores
  nps1 <- 100 * (sim1[,1] - sim1[,3])
  nps2 <- 100 * (sim2[,1] - sim2[,3])
  
  # Compute difference between nps scores
  diff.nps <- nps1 - nps2
  
  return(diff.nps)
}

# Calculate point estimate and credible intervals
nps.results.diff <- function(sim){
  # Point estimate and credible intervals
  point <- mean(sim)
  lb <- quantile(sim, 0.025)
  ub <- quantile(sim, 0.975)
  
  out <- c(lb, point, ub)
  return(out)
}

nps <- nps.sim.2test(9,0,2,2,0,9)
results <- nps.results(nps)
results

# Plot posterior
{hist(nps, breaks=100)
abline(v=results[[1]], col='red', lwd=3)
abline(v=mean(nps), col='red', lwd=3)
abline(v=results[[3]], col='red', lwd=3)}

################## POWER ANALYSIS #########################

# Length of one side of the credible intervals
# at a sample size of n
{
  ss <- c(20,30,50,75,100,150,200,300,500) # Sample sizes
  ci <- NULL # Credible intervals
  sim <- NULL # Results of simulation
  for(i in 1:length(ss)){
    sim <- nps.sim(ss[i]/2,ss[i]/2,ss[i])
    ci[i] <- (nps.results(sim)[[3]] - nps.results(sim)[[1]])/2
    print(ss[i])
  }
  pow <- data.frame(ss, ci)
  pow
}

# Plot of size of CIs at various sample sizes
{
plot(ss, ci, type='l', ylim=c(-25,25), ylab='Width of Confidence Interval of NPS',
     xlab='Sample Size', main='Power Analysis on NPS Confidence Interval')
lines(ss, -ci)
polygon(c(ss, rev(ss)), c(ci, rev(-ci)), col="skyblue", border = "red")
grid()
lines(c(ss[1],ss[length(ss)]), c(0,0) ,col='red')
for(i in 1:length(ss)) lines(c(ss[i],ss[i]),c(ci[i],-ci[i])) 
}

# Minimum difference needed to be detected 
# at a sample size of n (total sample size 2*n)
{
ss <- c(10,20,30,50,75,100,150,200,300,500) # Sample sizes
ci <- NULL
sim <- NULL
for(i in 1:length(ss)){
  sim <- nps.sim.2test(ss[i]/2,ss[i]/2,ss[i],ss[i]/2,ss[i]/2,ss[i])
  ci[i] <- (nps.results.diff(sim)[[3]] - nps.results.diff(sim)[[1]])/2
  print(i)
}
pow <- data.frame(ss, ci)
pow
}

# Plot of size of CIs of difference between two
# scores at various sample sizes
{
  plot(ss, ci, type='l', ylim=c(-51,51), ylab='Width of CI of Difference between Two NPS Scores',
       xlab='Sample Size', main='Power Analysis on Difference between Two NPS Scores')
  lines(ss, -ci)
  polygon(c(ss, rev(ss)), c(ci, rev(-ci)), col="lightgreen", border = "red")
  grid()
  lines(c(ss[1],ss[length(ss)]), c(0,0) ,col='red')
  for(i in 1:length(ss)) lines(c(ss[i],ss[i]),c(ci[i],-ci[i])) 
}

# Minimum difference needed to be detected 
# at a sample size of n (total sample size 2*n)
{
  ss1 <- c(10,20,30,50,75,100,150,200,300,500) # Sample size 1
  ss2 <- c(10,20,30,50,75,100,150,200,300,500) # Sample size 1
  ci <- data.frame()
  sim <- NULL
  for(i in 1:length(ss1)){
    for(j in 1:length(ss2)){
      sim <- nps.sim.2test(ss1[i]/2,ss1[i]/2,ss1[i],ss2[j]/2,ss2[j]/2,ss2[j])
      ci[i,j] <- (nps.results.diff(sim)[[3]] - nps.results.diff(sim)[[1]])/2
      print(c(i,j))
    }
  }
  row.names(ci) <- ss1
  names(ci) <- ss1
  ci
}

# Plot of size of CIs of difference between two
# scores at various sample sizes
{
  plot(ss, ci, type='l', ylim=c(-51,51), ylab='Width of CI of Difference between Two NPS Scores',
       xlab='Sample Size', main='Power Analysis on Difference between Two NPS Scores')
  lines(ss, -ci)
  polygon(c(ss, rev(ss)), c(ci, rev(-ci)), col="lightgreen", border = "red")
  grid()
  lines(c(ss[1],ss[length(ss)]), c(0,0) ,col='red')
  for(i in 1:length(ss)) lines(c(ss[i],ss[i]),c(ci[i],-ci[i])) 
}
