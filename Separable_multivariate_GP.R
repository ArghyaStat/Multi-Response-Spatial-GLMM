#See overleaf documentation for model

library(fields)
library(ggplot2)
library(viridis)
library(plot3D)
library(fBasics)  # For vectorize columns of a matrix
library(MCMCpack)
library(mvtnorm)

set.seed(3019)

# dimension of the random field
q <- 2

#as.integer(readline(prompt="Enter the dimension of the random field (q): "))

N <- 2e2 # No. of spatial locations

#simulating N locations over [0,1]^2 in locations matrix

locations <- matrix(0, nrow = N, ncol = 2)

locations[, 1] <- runif(N)
locations[, 2] <- runif(N)

# Adding a mean term

X = cbind(1, locations)

#design matrix in concatenated form of dim (Np * 1)
X_vec <- vec(t(X))

#Number of features in the spatial model
p <- ncol(X)

# True value of regression matrix beta
true_beta <- matrix(rnorm(p * q, mean = 0, sd = 2), nrow = p, ncol = q)

#beta_mat is the kronecker product of the coefficient matrix across locations

beta_mat <- diag(N) %x% t(true_beta)

# Fixed effect (I_n * B^T) X
mu_vec <- c(beta_mat %*% X_vec)

# True value of componentwise var-cov matrix
#true_Sigma <- matrix(runif(q^2, 0.5, 3), nrow = q) # Random symmetric positive definite matrix
#true_Sigma <- true_Sigma %*% t(true_Sigma)  # Ensuring positive definiteness

true_Sigma <- matrix(c(3,2,2,4), nrow = q, ncol =q, byrow = TRUE)

# true_phi

true_phi <- 0.1

# true_nu

true_nu <- 0.5

#true_r

true_r <- 0.8

# Calculates the Euclidean distance matrix

distmat <- rdist(locations)

# Calculates the correlation matrix from matern kernel

K <- Matern(distmat, range = true_phi, smoothness = true_nu)

# Calculates the separable covariance matrix with nugget effect

Omega <- (true_r*K  + (1-true_r)*diag(N)) %x% true_Sigma

# Generating the response vector of dimension Nq*1

Y_vec <- mu_vec + c(t(chol(Omega)) %*% rnorm(N*q))

# Saving necessary parameters and data
save(N,  p, q, locations, X, Y_vec, true_beta,
     true_Sigma, true_phi, true_nu, true_r, file = "separable_mgp_data.Rdata")

# MCMC Set up
# Log-likelihood function for Gaussian Process regression using Cholesky decomposition

log_likelihood <- function(beta, Sigma, phi, nu, r, Y_vec, locations) {
  
  N <- nrow(locations)
  
  X = cbind(1, locations)
  
  p <- ncol(X)
  
  #design matrix in concatenated form of dim (Np * 1)
  X_vec <- vec(t(X))
  
  distmat <- rdist(locations)
  
  if(phi > 0 & nu > 0 & r > 0 & r < 1){
    
    K.tilde <- r*Matern(distmat, range = phi, smoothness = nu) + (1- r)*diag(N)
    
    chol.K.tilde <- chol(K.tilde)
    prec.K.tilde <- chol2inv(chol.K.tilde)
    
    q <- ncol(Sigma)
    
    chol.Sigma <- chol(Sigma)
    prec.Sigma <- chol2inv(chol.Sigma)
    
    prec <- prec.K.tilde %x% prec.Sigma
    
    # Calculates the separable covariance matrix with nugget effect
    
    resid <- Y_vec - (diag(N) %x% t(beta)) %*% X_vec
    
    like.out <- - q*sum(log(diag(chol.K.tilde))) - N*sum(log(diag(chol.Sigma))) -
      0.5*sum((prec %*% resid) * resid)
    
  }else{
    
    like.out <- -Inf
  }
  
  return(like.out)
}

# log-posterior of the model parameters

log_posterior <- function(beta, Sigma, phi, nu, r, Y_vec, locations){
  
  prior_phi <- dunif(phi, 0, 1, log = TRUE)
  prior_nu <-  dlnorm(nu, meanlog = -1.2, sdlog = 1, log = TRUE)
  prior_Sigma <- log(diwish(Sigma, v =  q + 2, S = 100*diag(q)))
  prior_beta <- mvtnorm::dmvnorm(as.vector(vec(t(beta))), mean = rep(0, p*q),
                        sigma = 100*diag(p*q), log = TRUE)
  prior_r <- dunif(r, 0, 1, log = TRUE)
  
  return(log_likelihood(beta, Sigma, phi, nu, r, Y_vec, locations) +
           prior_beta +
           prior_Sigma +
           prior_phi +
           prior_nu +
           prior_r)
  
}

# Metropolis-Hastings algorithm for component-wise sampling

metropolis_hastings <- function(beta, Sigma, phi, nu, r, niters,
                                tuning_params) {
  
  beta_chain <- replicate(niters, matrix(0, p, q), simplify = F)
  Sigma_chain <- replicate(niters, matrix(0, q, q), simplify = F)
  phi_chain <- numeric(length = niters)
  nu_chain <- numeric(length = niters)
  r_chain <- numeric(length = niters)
  
  accept <- numeric(5)
  
  # Run Metropolis-Hastings
  for (iter in 2:niters) {
    
    if(iter %% ((niters)/10) == 0) print(paste0(100*(iter/(niters)), "%"))
    
    current_beta <- beta
    current_Sigma <- Sigma
    current_phi <- phi
    current_nu <- nu
    current_r <- r
    
    current_beta <- rnorm(p*q, beta, sqrt(tuning_params[1]))
    dim(current_beta) = c(p,q)
    
    # Compute log posterior for the proposed value
    log.r_beta <- log_posterior(current_beta,
                                Sigma,
                                phi,
                                nu,
                                r,
                                Y_vec, locations = locations) -
                  log_posterior(beta,
                                Sigma,
                                phi,
                                nu,
                                r,
                                Y_vec, locations = locations)
    
    # Accept or reject
    if (log(runif(1)) < log.r_beta) {
      
     
      beta <- current_beta
     
      accept[1] <- accept[1] + 1
      
    }
    
    current_Sigma <- riwish(v = q + 1 + tuning_params[2],
                            tuning_params[2]*Sigma)
    
    # Compute log posterior for the proposed value
    log.r_Sigma <- log_posterior(beta,
                                 current_Sigma,
                                 phi,
                                 nu,
                                 r,
                                 Y_vec, locations = locations) -
      log_posterior(beta,
                    Sigma,
                    phi,
                    nu,
                    r,
                    Y_vec, locations = locations) -
      log(diwish(W = current_Sigma,
                 v = q + 1 + tuning_params[2],
                 S = tuning_params[2]*Sigma)) +
      log(diwish(W = Sigma,
                 v = q + 1 + tuning_params[2],
                 S = tuning_params[2] * current_Sigma))
    
    # Accept or reject
    if (log(runif(1)) < log.r_Sigma) {
      
      Sigma <- current_Sigma
      
      accept[2] <- accept[2] + 1
      
    }
    
    current_phi <- rnorm(1, phi, sqrt(tuning_params[3]))
    
    # Compute log posterior for the proposed value
    log.r_phi <- log_posterior(beta,
                               Sigma,
                               current_phi,
                               nu,
                               r,
                               Y_vec, locations = locations) -
      log_posterior(beta,
                    Sigma,
                    phi,
                    nu,
                    r,
                    Y_vec, locations = locations)
    
    # Accept or reject
    if (log(runif(1)) < log.r_phi) {
      
      phi <- current_phi
      
      accept[3] <- accept[3] + 1
      
    }
    
    current_nu <- rnorm(1, nu, sqrt(tuning_params[4]))
    
    # Compute log posterior for the proposed value
    log.r_nu <- log_posterior(beta,
                              Sigma,
                              phi,
                              current_nu,
                              r,
                              Y_vec, locations = locations) -
      log_posterior(beta,
                    Sigma,
                    phi,
                    nu,
                    r,
                    Y_vec, locations = locations)
    
    # Accept or reject
    if (log(runif(1)) < log.r_nu) {
      
      nu <- current_nu
      
      accept[4] <- accept[4] + 1
      
    }
    
    current_r <- rnorm(1, r, sqrt(tuning_params[5]))
    #current_r <- runif(1, r - 0.5*tuning_params[5], sqrt(tuning_params[4]))
    
    # Compute log posterior for the proposed value
    log.r_r <- log_posterior(beta,
                             Sigma,
                             phi,
                             nu,
                             current_r,
                             Y_vec, locations = locations) -
               log_posterior(beta,
                             Sigma,
                             phi,
                             nu,
                             r,
                             Y_vec, locations = locations)
    
    # Accept or reject
    if (log(runif(1)) < log.r_r) {
      
      r <- current_r
      # accept overcounts
      accept[5] <- accept[5] + 1
      
    }
    
    beta_chain[[iter]] <- beta
    Sigma_chain[[iter]] <- Sigma
    phi_chain[iter] <- phi
    nu_chain[iter] <- nu
    r_chain[iter] <- r
    
  }
  
  theta_chain <- list("beta_sample" = beta_chain,
                      "Sigma_sample" = Sigma_chain,
                      "phi_sample" = phi_chain,
                      "nu_sample" = nu_chain,
                      "r_sample" = r_chain)
  
  # Calculate acceptance rate.
  print(paste("Acceptance = ", accept/niters))
  
  return(theta_chain)
  
}

# Sample initial parameters
beta <- matrix(rep(0, p*q), nrow = p, ncol = q)
Sigma <- diag(q)
phi <- 0.8
nu <- 1
r <- 0.5

# Number of iterations
niters <- 3e4

# Tuning parameters list

tuning_params <- c(0.02, 2e2, 5e-4, 1.5e-2, 1.5e-2)

# Run Metropolis-Hastings algorithm
theta_chain <- metropolis_hastings(beta, Sigma, phi, nu, r,
                                   niters = niters,
                                   tuning_params = tuning_params)

# Saving MCMC chain
save(theta_chain, file = "separable_mgp_MCMC_chain.Rdata")

# Traceplots

plot.ts(theta_chain$phi_sample, ylab = "phi", main = "Traceplot of phi")
abline(h = true_phi, col = 'blue', lwd = 2)

plot.ts(theta_chain$nu_sample, ylab = "nu", main = "Traceplot of nu")
abline(h = true_nu, col = 'blue', lwd = 2)

plot.ts(theta_chain$r_sample, ylab = "r", main = "Traceplot of r")
abline(h = true_r, col = 'blue', lwd = 2)

