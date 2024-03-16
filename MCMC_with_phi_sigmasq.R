library(fields)
library(plot3D)
library(geoR)

N <- 2e2 # No. of spatial locations 

locations <- matrix(0, nrow = N, ncol = 2)

locations[, 1] <- runif(N)
locations[, 2] <- runif(N) #simulating N locations


#plot(locations)
true_phi <- 0.1   #true value of phi
true_sigmasq <- 5


d <- rdist(locations) # Calculates the Euclidean distance matrix
cov <- true_sigmasq * exp(-d/true_phi) # Calculates the correlation matrix
Y <- c(t(chol(cov)) %*% rnorm(N)) # Draws a sample from multivariate normal

# MCMC Set up
# Log-likelihood function for Gaussian Process regression using Cholesky decomposition

log_likelihood <- function(theta_list, Y) {
  
  sigma_sq <- theta_list[[1]]
  phi <- theta_list[[2]]
  
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

# Target (the posterior of phi) with U(0,1) prior on phi

log_posterior_theta <- function(theta_list, Y){
  
  sigma_sq <- theta_list[[1]]
  phi <- theta_list[[2]]
  
  return(log_likelihood(theta_list, Y) + 
           dgamma(1/sigma_sq, shape = 0.01, rate = 0.01, log =TRUE) + 
           dunif(phi, 0, 1, log =TRUE))
}


# Metropolis-Hastings algorithm for component-wise sampling

metropolis_hastings <- function(theta_init, niters, proposal_sd) {
  
  n_params <- length(theta_init)
  theta_chain <- vector("list", length = n_params)
  error_list <- lapply(proposal_sd, function(h) rnorm(niters, 0, h))
  
  accept <- 0
  
  # Initialize chain
  for (j in 1:n_params) {
    
    theta_chain[[j]] <- numeric(niters)
    theta_chain[[j]][1] <- theta_init[[j]]
    
   # print(paste("Initialization Parameter:", j, "Initial value:", theta_chain[[j]][1]))
    
  }
 
  # Run Metropolis-Hastings
  for (iter in 2:niters) {
    
    if(iter %% ((niters)/10) == 0) print(paste0(100*(iter/(niters)), "%"))
    
    for (j in 1:n_params) {
      
      # At this level current takes the theta_chain at (iter-1)
      current <- lapply(theta_chain, function(x) x[iter - 1])
    
      
      # Propose a new value for parameter j
      proposed_theta <- theta_chain[[j]][iter - 1] + error_list[[j]][iter]
      
      
      # Just update 'only" the jth component of theta_chain
      
      current[[j]] <- proposed_theta
      
      # Compute log posterior for the proposed value
      
      
      log.r <- log_posterior_theta(current, Y) - 
               log_posterior_theta(lapply(theta_chain, function(x) x[iter - 1]), 
                                   Y)
      
      # Print diagnostic information
     # print(paste("Iteration:", iter, "Parameter:", j, "Proposed theta:", proposed_theta, "Log-ratio:", log.r))
      
      # Accept or reject
      if (log(runif(1)) < log.r) {
        
        
      # suspected wrong step  
       theta_chain[[j]][iter] <- current[[j]]
       
       # accept overcounts
       accept <- accept + 1
        
      } else {
        
       theta_chain[[j]][iter] <- theta_chain[[j]][iter - 1]
        
      }
    }
    # Print theta_chain values for each iteration
    #print(paste("Theta_chain at iteration", iter, ":", theta_chain))
  }
  
  # Calculate acceptance rate which is > 1 wrong.
  print(paste("Acceptance = ", accept/niters))
  
  return(theta_chain)
}


# Sample initial parameters
theta_init <- list(sigma_sq = 12, phi = 0.5)  # Ensure the correct order and naming

# Number of iterations
niters <- 10000

# Standard deviation for proposal distribution
proposal_sd <- c(0.2, 0.05)

# Run Metropolis-Hastings algorithm
theta_chain <- metropolis_hastings(theta_init = theta_init, niters = niters, proposal_sd = proposal_sd)
