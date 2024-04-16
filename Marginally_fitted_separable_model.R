# Marginal model is fitted with Sigma = I_q

library(fields)
library(ggplot2)
library(viridis)
library(plot3D)
library(fBasics)  # For vectorize columns of a matrix
library(MCMCpack)
library(mvtnorm)
library(SimTools)

set.seed(3019)

# dimension of the random field
q <- 2

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

log_likelihood <- function(beta, phi, nu, r, Y_vec, locations) {
  
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
    
    q <- ncol(beta)
    
    
    prec <- prec.K.tilde %x% diag(q)
    
    # Calculates the separable covariance matrix with nugget effect
    
    resid <- Y_vec - (diag(N) %x% t(beta)) %*% X_vec
    
    like.out <- - q*sum(log(diag(chol.K.tilde))) - 0.5*sum((prec %*% resid) * resid)
    
  }else{
    
    like.out <- -Inf
  }
  
  return(like.out)
}

# log-posterior of the model parameters

log_posterior <- function(beta, phi, nu, r, Y_vec, locations){
  
  prior_phi <- dunif(phi, 0, 1, log = TRUE)
  prior_nu <-  dlnorm(nu, meanlog = -1.2, sdlog = 1, log = TRUE)
  prior_beta <- mvtnorm::dmvnorm(as.vector(vec(t(beta))), mean = rep(0, p*q),
                                 sigma = 100*diag(p*q), log = TRUE)
  prior_r <- dunif(r, 0, 1, log = TRUE)
  
  return(log_likelihood(beta, phi, nu, r, Y_vec, locations) +
           prior_beta +
           prior_phi +
           prior_nu +
           prior_r)
  
}

# Metropolis-Hastings algorithm for component-wise sampling

metropolis_hastings <- function(beta, phi, nu, r, niters,
                                tuning_params) {
  
  beta_chain <- replicate(niters, matrix(0, p, q), simplify = F)
  phi_chain <- numeric(length = niters)
  nu_chain <- numeric(length = niters)
  r_chain <- numeric(length = niters)
  
  accept <- numeric(4)
  
  # Run Metropolis-Hastings
  for (iter in 2:niters) {
    
    if(iter %% ((niters)/10) == 0) print(paste0(100*(iter/(niters)), "%"))
    
    current_beta <- beta
    current_phi <- phi
    current_nu <- nu
    current_r <- r
    
    current_beta <- rnorm(p*q, beta, sqrt(tuning_params[1]))
    dim(current_beta) = c(p,q)
    
    # Compute log posterior for the proposed value
    log.r_beta <- log_posterior(current_beta,
                                phi,
                                nu,
                                r,
                                Y_vec, locations = locations) -
      log_posterior(beta,
                    phi,
                    nu,
                    r,
                    Y_vec, locations = locations)
    
    # Accept or reject
    if (log(runif(1)) < log.r_beta) {
      
      
      beta <- current_beta
      
      accept[1] <- accept[1] + 1
      
    }
    
    
    current_phi <- rnorm(1, phi, sqrt(tuning_params[2]))
    
    # Compute log posterior for the proposed value
    log.r_phi <- log_posterior(beta,
                               current_phi,
                               nu,
                               r,
                               Y_vec, locations = locations) -
      log_posterior(beta,
                    phi,
                    nu,
                    r,
                    Y_vec, locations = locations)
    
    # Accept or reject
    if (log(runif(1)) < log.r_phi) {
      
      phi <- current_phi
      
      accept[2] <- accept[2] + 1
      
    }
    
    current_nu <- rnorm(1, nu, sqrt(tuning_params[3]))
    
    # Compute log posterior for the proposed value
    log.r_nu <- log_posterior(beta,
                              phi,
                              current_nu,
                              r,
                              Y_vec, locations = locations) -
      log_posterior(beta,
                    phi,
                    nu,
                    r,
                    Y_vec, locations = locations)
    
    # Accept or reject
    if (log(runif(1)) < log.r_nu) {
      
      nu <- current_nu
      
      accept[3] <- accept[3] + 1
      
    }
    
    current_r <- rnorm(1, r, sqrt(tuning_params[4]))
    #current_r <- runif(1, r - 0.5*tuning_params[5], sqrt(tuning_params[4]))
    
    # Compute log posterior for the proposed value
    log.r_r <- log_posterior(beta,
                             phi,
                             nu,
                             current_r,
                             Y_vec, locations = locations) -
      log_posterior(beta,
                    phi,
                    nu,
                    r,
                    Y_vec, locations = locations)
    
    # Accept or reject
    if (log(runif(1)) < log.r_r) {
      
      r <- current_r
      # accept overcounts
      accept[4] <- accept[4] + 1
      
    }
    
    beta_chain[[iter]] <- beta
    phi_chain[iter] <- phi
    nu_chain[iter] <- nu
    r_chain[iter] <- r
    
  }
  
  theta_chain <- list("beta_sample" = beta_chain,
                      "phi_sample" = phi_chain,
                      "nu_sample" = nu_chain,
                      "r_sample" = r_chain)
  
  # Calculate acceptance rate.
  print(paste("Acceptance = ", accept/niters))
  
  return(theta_chain)
  
}

# Sample initial parameters
beta <- matrix(rep(0, p*q), nrow = p, ncol = q)
phi <- 0.5
nu <- 1
r <- 0.2

# Number of iterations
niters <- 2e4

# Tuning parameters list

tuning_params <- c(0.08, 1e-3, 1.5e-2, 1e-2)

# Run Metropolis-Hastings algorithm
theta_chain <- metropolis_hastings(beta, phi, nu, r,
                                   niters = niters,
                                   tuning_params = tuning_params)

# Saving MCMC chain
save(theta_chain, file = "wrong_separable_mgp_MCMC_chain.Rdata")

theta_chain <- load(file = "separable_mgp_MCMC_chain.Rdata")

# Traceplots

par(mfrow = c(1,1))
trace_phi <- plot.ts(theta_chain$phi_sample, ylab = "phi", main = "Traceplot of phi")
abline(h = true_phi, col = 'blue', lwd = 2)

trace_nu <- plot.ts(theta_chain$nu_sample, ylab = "nu", main = "Traceplot of nu")
abline(h = true_nu, col = 'blue', lwd = 2)

trace_r <- plot.ts(theta_chain$r_sample, ylab = "r", main = "Traceplot of r")
abline(h = true_r, col = 'blue', lwd = 2)

# acfplots

acf_phi <- acf(theta_chain$phi_sample, main = "ACF plot of phi", lag.max = 200)
acf_nu <- acf(theta_chain$nu_sample, main = "ACF plot of nu", lag.max = 200)
acf_r <- acf(theta_chain$r_sample, main = "ACF plot of r", lag.max = 200)


## Output analysis using SimTools
traceplot(theta_chain$phi_sample, ylab = "phi",
          main = "Trace plot of phi")
abline(h = true_phi, col = "blue")
acfplot(theta_chain$phi_sample, lag.max = 200, main = "ACF plot of phi")
out_phi <- as.Smcmc(theta_chain$phi_sample)
plot(out_phi)


out_nu <- as.Smcmc(theta_chain$nu_sample)
plot(out_nu)

out_r <- as.Smcmc(theta_chain$r_sample)
plot(out_r)




labels_beta <- list()
for (i in 1:p) {
  for (j in 1:q) {
    labels_beta <- c(labels_beta, list(paste0('beta (', i, ',', j, ')')))
  }
}



par(mfrow=c(p,q))

for (i in 1:p) {
  for (j in 1:q) {
    plot(1:niters, sapply(theta_chain$beta_sample, function(x) x[i, j]), 
         type='l', col=1, xlab='Iteration', ylab='beta',
         main=labels_beta[(i-1)*q + j])
    abline(h = true_beta[i,j], col = 'blue', lwd = 2)
    
  }
}

par(mfrow=c(p,q))

for (i in 1:p) {
  for (j in 1:q) {
    
    acf(sapply(theta_chain$beta_sample, function(x) x[i, j]), 
        main = labels_beta[(i-1)*q + j], lag.max = 200)
    
  }
}


# Model Comparison

DIC <- function(theta, niters, burnIn, Y_vec, locations) {
  
  
  
  log_like_samples <- numeric(niters - burnIn)
  
  post_mean_beta <- matrix(0, p, q)
  post_mean_phi <- numeric(1)
  post_mean_nu <- numeric(1)
  post_mean_r <- numeric(1)
  
  for(iter in (burnIn + 1):niters){
    
    log_like_samples[iter] <- log_likelihood(theta$beta_sample[[iter]],
                                             theta$phi_sample[[iter]],
                                             theta$nu_sample[[iter]],
                                             theta$r_sample[[iter]],
                                             Y_vec,
                                             locations)
    
    
    post_mean_beta <- post_mean_beta + theta$beta_sample[[iter]]
    post_mean_phi <- post_mean_phi + theta$phi_sample[[iter]]
    post_mean_nu <- post_mean_nu + theta$nu_sample[[iter]]
    post_mean_r <- post_mean_r + theta$r_sample[[iter]]
    
  }
  
  post_mean_beta <- post_mean_beta/(niters - burnIn)
  post_mean_phi <- post_mean_phi/(niters - burnIn)
  post_mean_nu <- post_mean_nu/(niters - burnIn)
  post_mean_r <- post_mean_r/(niters - burnIn)
  
  # Compute the mean of the log-likelihood samples
  log_like_mean <- mean(log_like_samples)
  
  # likelihood at posterior mean
  log_like_post <- log_likelihood(post_mean_beta,
                                  post_mean_phi,
                                  post_mean_nu,
                                  post_mean_r,
                                  Y_vec = Y_vec,
                                  locations = locations)
  
  # Compute the Deviance
  DIC <- -4 * log_like_mean + 2* log_like_post
  
  return(DIC)
}


DIC_marginal <- DIC(theta = theta_chain, 
                 niters = 2e4, 
                 burnIn = 1e4,
                 Y_vec = Y_vec,
                 locations = locations)

