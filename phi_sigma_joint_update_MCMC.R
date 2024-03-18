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



log_posterior_theta <- function(theta_vec, Y){
  
  sigma_sq <- theta_vec[1]
  phi <- theta_vec[2]

  # Target with U(0,1) prior on phi and IG(0.01, 0.01) on sigma_sq
  
  prior_sigma_sq <- dgamma(1/sigma_sq, shape = 0.01, rate = 0.01, log =TRUE)
  prior_phi <- dunif(phi, 0, 1, log =TRUE))
  
  return(log_likelihood(theta_vec, Y) + 
           prior_sigma_sq + 
           prior_phi
}


# Metropolis-Hastings algorithm for joint sampling
metropolis_hastings <- function(theta_init, niters, proposal_sd) {
  
  n_params <- length(theta_init)
  theta_chain <- array(dim = c(niters, n_params))
  
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
niters <- 1e4

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
         
         
# Joint contour plot

library(ggplot2)
library(viridis)

# Create a grid of sigma_sq and phi values
sigma_sq_values <- seq(0.001, 3, length.out = 100)  # Adjust the range and resolution as needed
phi_values <- seq(0, 0.25, length.out = 100)           # Adjust the range and resolution as needed
grid <- expand.grid(sigma_sq = sigma_sq_values, phi = phi_values)

# Calculate log posterior density for each combination of sigma_sq and phi
log_posterior_values <- apply(grid, 1, function(params) {
  log_posterior_theta(params, Y)
})

# Combine grid with log posterior values
grid$log_posterior <- log_posterior_values

# Convert log posterior values to actual density values
grid$density <- exp(grid$log_posterior)

# Create contour plot
p <- ggplot(grid, aes(sigma_sq, phi, z = density)) +
  geom_contour_filled(aes(fill = ..level..), bins = 20) +
  labs(x = expression(sigma^2), y = expression(phi)) +
  ggtitle("Joint Posterior Density of sigma_sq and phi") +
  theme_minimal() +
  theme(legend.position = "right")

# Show the plot
p
