
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
  
  # Initialize chain
  for (j in 1:n_params) {
    
    theta_chain[[j]] <- numeric(niters)
    theta_chain[[j]][1] <- theta_init[[j]]
    
    print(paste("Initialization Parameter:", j, "Initial value:", theta_chain[[j]][1]))
    
  }
  
  accept <- 0
  
  # Run Metropolis-Hastings
  for (iter in 2:niters) {
    
    #if(i %% ((iters)/10) == 0) print(paste0(100*(i/(iters)), "%"))
    
    for (j in 1:n_params) {
      
      # Propose a new value for parameter j
      proposed_theta <- theta_chain[[j]][iter - 1] + error_list[[j]][iter]
      
      # Compute log posterior for the proposed value
      
      
      log.r <- log_posterior_theta(proposed_theta, Y) - 
               log_posterior_theta(theta_chain[[j]][iter - 1], Y)
      
      # Print diagnostic information
      print(paste("Iteration:", iter, "Parameter:", j, "Proposed theta:", proposed_theta, "Log-ratio:", log.r))
      
      # Accept or reject
      if (log(runif(1)) < log.r) {
        
        theta_chain[[j]][iter] <- proposed_theta
        accept <- accept + 1
        
      } else {
        
        theta_chain[[j]][iter] <- theta_chain[[j]][iter - 1]
        
      }
    }
    # Print theta_chain values for each iteration
    print(paste("Theta_chain at iteration", iter, ":", theta_chain))
  }
  
  # Calculate acceptance rate
  print(paste("Acceptance = ", accept/niters))
  
  return(theta_chain)
}


# Sample initial parameters
theta_init <- list(sigma_sq = 20, phi = 0.8)  # Ensure the correct order and naming

# Number of iterations
niters <- 2

# Standard deviation for proposal distribution
proposal_sd <- c(0.2, 0.05)

# Run Metropolis-Hastings algorithm
theta_chain <- metropolis_hastings(theta_init = theta_init, niters = niters, proposal_sd = proposal_sd)

# Run Metropolis-Hastings algorithm
theta_chain <- metropolis_hastings(theta_init = theta_init, niters = niters, proposal_sd = proposal_sd)


# Output Analysis

# Calculate the auto-correlation function (ACF) for each parameter
calculate_acf <- function(samples) {
  acf_values <- array(0, dim = c(length(samples), length(samples[[1]])))
  for (i in 1:length(samples)) {
    for (j in 1:length(samples[[1]])) {
      acf_values[i, j] <- acf(samples[[i]][[j]], plot = FALSE)$acf[2]
    }
  }
  return(acf_values)
}

# Effective sample size (ESS) calculation for each parameter
calculate_ess <- function(acf_values) {
  n <- length(acf_values)
  ess_values <- array(0, dim = c(length(acf_values[[1]])))
  for (i in 1:length(acf_values[[1]])) {
    rho <- 1 + 2 * sum(acf_values[, i])
    ess_values[i] <- n / rho
  }
  return(ess_values)
}


# Function to create density plots for each component of the MCMC chain
density_plots <- function(theta_chain, main = "Density Plots") {
  par(mfrow = c(length(theta_chain), 1))
  for (i in 1:length(theta_chain)) {
    density <- density(theta_chain[[i]])
    plot(density, main = paste("Component", i), xlab = "", ylab = "")
  }
  par(mfrow = c(1, 1))  # Reset plotting layout
  title(main = main)
}

# Function to create ACF plots for each component of the MCMC chain
acf_plots <- function(theta_chain, main = "ACF Plots") {
  par(mfrow = c(length(theta_chain), 1))
  for (i in 1:length(theta_chain)) {
    acf_theta <- acf(theta_chain[[i]], main = paste("Component", i))
  }
  par(mfrow = c(1, 1))  # Reset plotting layout
  title(main = main)
}

# Function to create trace plots for each component of the MCMC chain
trace_plots <- function(theta_chain, main = "Trace Plots") {
  par(mfrow = c(length(theta_chain), 1))
  for (i in 1:length(theta_chain)) {
    plot(theta_chain[[i]], type = "l", xlab = "Iteration", ylab = paste("Component", i))
  }
  par(mfrow = c(1, 1))  # Reset plotting layout
  title(main = main)
}

# Function to compute effective sample size (ESS) for each component of the MCMC chain
effective_sample_size <- function(theta_chain) {
  n_params <- length(theta_chain)
  ess <- numeric(n_params)
  for (i in 1:n_params) {
    acf_theta <- acf(theta_chain[[i]], plot = FALSE)
    rho_k <- acf_theta$acf[-1]  # Exclude lag 0 autocorrelation
    n <- length(rho_k)
    var_estimate <- 1 + 2 * sum(rho_k)
    ess[i] <- n / var_estimate
  }
  return(ess)
}


# Assuming theta_chain is the output of the MCMC algorithm
# Call the functions to create plots
density_plots(theta_chain, main = "Density Plots")
acf_plots(theta_chain, main = "ACF Plots")
trace_plots(theta_chain, main = "Trace Plots")
ess_values <- effective_sample_size(theta_chain)
print(ess_values)