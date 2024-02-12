rm(list = ls())

set.seed(3019)


mydir <- "C:/Users/Arghya/Documents/GitHub/Multi-Response-Spatial-GLMM"
setwd(mydir)

# Simulating with  Phi = 0.05

library(fields)
library(plot3D)
library(geoR)

N <- 2e2 # No. of spatial locations

loc <- cbind(runif(N), runif(N)) #simulating 100 locations
#plot(loc)
true_phi <- 0.1  #rep(0.05,100)  #true value of phi


d <- rdist(loc) # Calculates the Euclidean distance matrix
cov <- exp(-d/true_phi) # Calculates the correlation matrix
Y <- c(t(chol(cov)) %*% rnorm(N)) # Draws a sample from multivariate normal

scatter2D(x = loc[ , 1], y = loc[ , 2], colvar = Y, pch = 15)


# MCMC Set up

loglike_phi <- function(phi, Y){
    
    #nu <- 0.5
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

# Target

target_phi <- function(phi, Y){
  
  return(loglike_phi(phi, Y) + dunif(phi, 0, 1, log = T))
}

# Metropolis-Hastings Algorithm

mh_phi <- function(phi.init, iters, burnin){
  
  phi.chain <- array(dim = c(iters + burnin, 1))
  phi.chain[1] <- phi.init 
  
  for (i in 2:(iters + burnin))
  {
    
    proposal <- rnorm(n = 1, mean = phi.chain[i-1], sd = 0.1) 
    
    r = exp(target_phi(proposal, Y) - target_phi(phi.chain[i-1], Y))
    
    
    acc_ratio = min(r,1)
    
    
    if(runif(1) < acc_ratio){
      
      phi.chain[i] = proposal
      
    } else {
      
      phi.chain[i] = phi.chain[i-1]
      
    }
    
  }
  
  acceptance = 1 - mean(duplicated(phi.chain[-(1:burnin)]))
  
  return(phi.chain)
  
}

MCMC.out <- mh_phi(phi.init = 0.5, iters <- 2e3, burnin <- 0)

save(MCMC.out, file = "fit.observed.GP.Rdata")


# Plotting the MCMC output
par(mfrow = c(1, 2))

pdf("Output_plots.pdf", width = 4, height = 4)

p1 <- plot(MCMC.out, xlab = "Iteration", ylab = "phi", type = "l", 
          main = "Trace plot of phi")
abline(h = true_phi, col = 'blue', lwd = 2)

p2 <- acf(MCMC.out)



#---------------------



plot(fit.observed.GP$mu, type = "l")
plot(fit.observed.GP$epsilon.bar, type = "l")

acf(fit.observed.GP$mu)
acf(fit.observed.GP$epsilon.bar)

dev.off()



