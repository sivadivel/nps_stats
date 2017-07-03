# NPS Statistical Analysis
# Levi Davis 6/15/2017

# Housekeeping
{
  # Turn off scientific notation
  options(scipen=999)
  }

#################### STATS ###########################

# Function giving NPS with 95% CIs
# np is the number of promoters
# nd is the number of detractors
# n is the total number of respondents
nps.ci <- function(np,nd,n) {
  
  p <- np/n  # Proportion of promoters
  d <- nd/n  # Proportion of detractors
  
  nps <- 100 * (p-d) # NPS
  
  # Standard Error of Net Promoter Score
  # Source: http://www.cybaea.net/journal/2016/01/14/NPS-confidence-intervals-stats/
  vr <- (p*(1-p) + d*(1-d) + 2*p*d) / n
  se <- 100 * sqrt(vr)
  
  # Lower and upper bound of NPS
  lb <- nps - 1.96 * se 
  ub <- nps + 1.96 * se
  
  # Keeping bounds between -100 and 100
  for(i in 1:length(lb)) if (lb[i] < -100) lb[i] <- -100
  for(i in 1:length(ub)) if (ub[i] > 100) ub[i] <- 100
  
  # Rounding to two decimal places
  nps <- round(nps, 2)
  lb <- round(lb, 2)
  ub <- round(ub, 2)
  
  # Results
  out <- list(lb, nps, ub)
  names(out) <- c('Lower Bound', 'Point Estimate', 'Upper Bound')
  return(out)
}

nps.ci(9,2,11)

# Function testing difference between two sets of respodents
# Returns lower bound, point estimate and upper bound 
# and the p-value of difference between the two NPS scores.
# np1 and np2 are the number of promoters in group one and two.
# nd1 and nd2 are the number of detractors in group one and two.
# n1 and n2 are the number of total respondents in group one and two.
nps.test <- function(np1,nd1,n1,np2,nd2,n2) {
  
  p1 <- np1/n1  # Proportion of promoters group 1
  d1 <- nd1/n1  # Proportion of detractors group 1
  p2 <- np2/n2  # Proportion of promoters group 2
  d2 <- nd2/n2  # Proportion of detractors group 2
  
  nps1 <- 100 * (p1-d1) # NPS number 1
  nps2 <- 100 * (p2-d2) # NPS number 2
  
  
  d <- nps1 - nps2 # Difference between group 1 and 2 scores
  
  # Standard Error of difference between group 1 and 2 scores
  se <- sqrt(((p1*(1-p1) + d1*(1-d1) + 2*p1*d1)*100^2 / n1) +
              (p2*(1-p2) + d2*(1-d2) + 2*p2*d2)*100^2 / n2)
  
  # Lower and upper bound of difference between two groups
  lb <- d - 1.96*se
  ub <- d + 1.96*se
  
  # Round to zero decimal places
  lb <- round(lb, 0)
  ub <- round(ub, 0)
  d <- round(d, 0)
  
  # Calculate p-value
  pval <- 2 * (1 - pnorm(abs(d), sd=se))
  
  # Errors
  if (np1+nd1>n1 || np2+nd2>n2) return("Promoters + Detractors > Total Respondents")
  
  # Results
  out <- list(c(lb, d, ub), pval,d)
  names(out) <- c('Point estimate with CIs', 'p-value', 'Difference')
  return(list(c(lb, d, ub), pval,d))
}

nps.test(20,10,30,12,18,30)

################# POWER ANALYSIS #####################

# Width of one side of confidence intervals
# Finds width of one side CIs given NPS between n1 and n2
nps.power.ci <- function(p, d, n1, n2){
  s <- n1:n2  # Vector of values counting up from smallest to largest 
  np <- p * s # Number of promoters
  nd <- d * s # Number of detractors
  
  # Diference between the upper and lower bounds
  wci <- 0.5 * (nps.ci(np, nd, s)[[3]] - nps.ci(np, nd, s)[[1]])
  return(wci)
}

# Plot width of confidence interval
{
plot(nps.power.ci(0.5, 0.5, 10, 100), type='l', ylim=c(0,70), xlab='Number of Responents',
     ylab='Width of one side of confidence interval')
lines(nps.power.ci(0.75, 0.25, 10, 100), col='red')
lines(nps.power.ci(0.25, 0.25, 10, 100), col='blue')
grid()
}

N <- 1000
# (np1,nd1,n1,np2,nd2,n2)
nps.power.test2 <- function(N){
  d <- (1:100)/100
  n <- N
  np1 <- n/2
  nd1 <- n/2
  np2 <- n/2 * (1+d)
  nd2 <- n/2 * (1-d)
  test <- nps.test(np1,nd1,n,np2,nd2,n)
  pvals <- test[[2]]
  dif <- abs(test[[3]])
  return(list(pvals, dif))
}

plot(nps.power.test2(10)[[2]], nps.power.test2(10)[[1]], type='l', xlab='Difference',
     ylab='P-value', ylim=c(0,0.25))
rect(0,0,100,0.05, col='skyblue')
lines(nps.power.test2(10)[[2]], nps.power.test2(10)[[1]])
lines(nps.power.test2(20)[[2]], nps.power.test2(20)[[1]], col='red')
lines(nps.power.test2(30)[[2]], nps.power.test2(30)[[1]], col='blue')

nps.power.test2(10)

d <- (1:100)/100
n <- 30
np1 <- n/2
nd1 <- n/2
np2 <- n/2 * (1+d)
nd2 <- n/2 * (1-d)
np2 - nd2