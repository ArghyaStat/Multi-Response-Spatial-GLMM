library(fields)
library(plot3D)
library(geoR)


N <- 2e2 # No. of spatial locations 

locations <- matrix(0, nrow = N, ncol = 2)

locations[, 1] <- runif(N)
locations[, 2] <- runif(N) #simulating N locations


#plot(locations)
true_phi <- 0.1  #true value of phi
true_sigmasq <- 2  #true value of sigma_sq


d <- rdist(locations) # Calculates the Euclidean distance matrix
cov <- true_sigmasq * exp(-d/true_phi) # Calculates the correlation matrix
Y <- c(t(chol(cov)) %*% rnorm(N)) # Draws a sample from multivariate normal


# MCMC Set up
# Log-likelihood function for Gaussian Process regression using Cholesky decomposition

log_likelihood <- function(theta_vec, Y) {
  
  sigma_sq <- theta_vec[1]
  phi <- theta_vec[2]
  
  n <- length(Y)
  
  if(phi > 0 && sigma_sq > 0){
    
    cov_mat <- sigma_sq * exp(-d/phi)
    chol.cov <- chol(cov_mat)
    prec <- chol2inv(chol.cov)
    like.out <-  -sum(log(diag(chol.cov))) - 0.5*sum((prec %*% Y) * Y)
  }else{
    
    like.out <- -Inf
  }
  
  return(like.out)
}

# Target with U(0,1) prior on phi and IG(0.01, 0.01) on sigma_sq

log_posterior_theta <- function(theta_vec, Y){
  
  sigma_sq <- theta_vec[1]
  phi <- theta_vec[2]
  
  return(log_likelihood(theta_vec, Y) + 
           dgamma(1/sigma_sq, shape = 0.01, rate = 0.01, log =TRUE) + 
           dunif(phi, 0, 1, log =TRUE))
}


# Metropolis-Hastings algorithm for joint sampling
metropolis_hastings <- function(theta_init, niters, proposal_sd) {
  
  n_params <- length(theta_init)
  theta_chain <- array(dim = c(niters, 2))
  
  error_mat <- cbind(rnorm(niters, 0, proposal_sd[1]), rnorm(niters, 0, proposal_sd[2]))
  accept <- 0
  
  theta_chain[1,] <- theta_init
  
  for (i in 2:niters) {
    
    
    if(i %% ((niters)/10) == 0) print(paste0(100*(i/(niters)), "%"))
    
    
    # Propose new value for the parameter
    proposal_theta <- theta_chain[i-1,] + error_mat[i,]
    
    log.r <- log_posterior_theta(proposal_theta, Y) - 
      log_posterior_theta(theta_chain[i-1,] , Y)
    
    # Accept or reject
    if (log(runif(1)) < log.r) {
      
      
      # suspected wrong step  
      theta_chain[i,] <- proposal_theta
      
      # accept counter
      accept <- accept + 1
      
    } else {
      
      theta_chain[i,] <- theta_chain[i-1,]
      
    }
  }
  
  # Calculate acceptance rate 
  print(paste("Acceptance = ", accept/niters))
  
  return(theta_chain)
}



# Sample initial parameters
theta_init <- c(10,0.5)  # Ensure the correct order and naming

# Number of iterations
niters <- 10000

# Standard deviation vector for proposal distribution
proposal_sd <- c(0.5, 0.08)

# Run Metropolis-Hastings algorithm
theta_chain <- metropolis_hastings(theta_init = theta_init, niters = niters, proposal_sd = proposal_sd)


# Traceplots

par(mfrow = c(ncol(theta_chain), 1))


plot.ts(theta_chain[,1], ylab = "sigma_sq", main = "Traceplot of sigma_sq")
abline(h = true_sigmasq, col = 'blue', lwd = 2)

plot.ts(theta_chain[,2], ylab = "phi", main = "Traceplot of phi")
abline(h = true_phi, col = 'blue', lwd = 2)