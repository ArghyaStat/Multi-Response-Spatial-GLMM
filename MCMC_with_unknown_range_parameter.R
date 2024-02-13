###################################
## No covariate model, with only
## unknown range parameter
## Look at ----- overleaf for theory
###################################

set.seed(3019)



# Simulating with  Phi = 0.1

library(fields)
library(plot3D)
library(geoR)
library(SimTools)

N <- 2e2 # No. of spatial locations 

loc <- matrix(0, nrow = N, ncol = 2)
loc[, 1] <- runif(N)
loc[, 2] <- runif(N) #simulating 100 locations
plot(loc)


true_phi <- 0.1  #true value of phi


d <- rdist(loc) # Calculates the Euclidean distance matrix
cov <- exp(-d/true_phi) # Calculates the correlation matrix
Y <- c(t(chol(cov)) %*% rnorm(N)) # Draws a sample from multivariate normal

# plotting bivariate density on the domain [0,1]^2

scatter2D(x = loc[ , 1], y = loc[ , 2], colvar = Y, pch = 15) 


# MCMC Set up


# log likelihood of phi using cholesky decomposition

loglike_phi <- function(phi, Y){
    
    if(phi > 0){
      
      cov_mat <- exp(-d/phi)
      chol.cov <- chol(cov_mat)
      prec <- chol2inv(chol.cov)
      like.out <-  -sum(log(diag(chol.cov))) - 0.5*sum((prec %*% Y) * Y)
    }else{
       
      like.out <- -Inf
    }
    
    return(like.out)

}

# Target (the posterior of phi) with U(0,1) prior on phi

target_phi <- function(phi, Y){
  
  return(loglike_phi(phi, Y) + dunif(phi, 0, 1, log = TRUE))
}

# Metropolis-Hastings Algorithm

mh_phi <- function(phi.init, iters, burnin = 0, h = 0.1){
  
  phi.chain <- numeric(length = (iters + burnin) )
  phi.chain[1] <- phi.init 
  accept <- 0
  errors <- rnorm(iters + burnin, sd = h)
  for (i in 2:(iters + burnin))
  {
    
    # Gaussian proposal with SD = 0.1
    if(i %% ((iters+burnin)/10) == 0) print( paste0(100*(i/(iters + burnin)), "%"))
    
    proposal <- phi.chain[i-1] + errors[i]
    log.r = target_phi(proposal, Y) - target_phi(phi.chain[i-1], Y) 
    
    if(log(runif(1)) < log.r){
      
      accept <- accept + 1
      phi.chain[i] = proposal
      
    } else {
      
      phi.chain[i] = phi.chain[i-1]
      
    }
  }
  print(paste("Acceptance = ", (accept/(iters + burnin))))
  return(phi.chain)
}

MCMC.out <- mh_phi(phi.init = 0.1, iters = 2e3, h = .1)

# acceptance rate

(acceptance = 1 - mean(duplicated(MCMC.out)))


# Plotting the MCMC output

par(mfrow = c(1, 2))

p1 <- plot(MCMC.out, xlab = "Iteration", ylab = "phi", type = "l", 
          main = "Trace plot of phi")
abline(h = true_phi, col = 'blue', lwd = 2)

p2 <- acf(MCMC.out)

## Output analysis using SimTools
traceplot(MCMC.out)
abline(h = true_phi, col = "blue")
acfplot(MCMC.out)
out <- as.Smcmc(MCMC.out)
plot(out)

#---------------------
