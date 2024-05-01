##### Metropolis-Gibbs Algorithm for updating params in 2 dimensional continuous-response separable GP model
##### Parameters in the model: B_{p*q}, Sigma_{q*q}, phi, nu, r
##### Full conditional distributions of B and Sigma are documented in overleaf section 3.3

ibrary(fields)
library(rlist)
library(mniw)
library(fBasics)  # For vectorize columns of a matrix
library(MCMCpack)
library(mvtnorm)
library(SimTools)

# For visualization
library(ggplot2)
library(viridis)
library(plot3D)


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
#save(N,  p, q, locations, X, Y_vec, true_beta,
     # true_Sigma, true_phi, true_nu, true_r, file = "separable_mgp_data.Rdata")

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


# Define function for matrix-normal density
dmatnorm <- function(X, M, U, V, log = FALSE) {
  
  p <- nrow(X)
  q <- ncol(X)
  
  denom <- -0.5 * p * q * log(2 * pi)  
  -0.5 * (p * sum(log(diag(chol(V)))) + q * sum(diag(log(chol(U)))))
  
  if (!isSymmetric(V)) {
    stop("V must be symmetric positive definite.")
  }
  if (!isSymmetric(U)) {
    stop("U must be symmetric positive definite.")
  }
  
  V.prec <- chol2inv(chol(V))
  U.prec <- chol2inv(chol(U))
  
  exponent <- -0.5 * sum(V.prec %*% diag(t(X - M) %*% U.prec %*% (X-M)))
  
  if (log) {
    return(denom + exponent)
  } else {
    return(exp(denom + exponent))
  }
}


# Define function for inverse-Wishart density

dinvwish <- function(A, df, S, log = FALSE) {
  
  q <- nrow(A)
  
  const <- (0.5 * q * df) * sum(log(diag(chol(S)))) 
  - (0.5 * df) * log(2) - lgamma(0.5 * q * df ) - 0.5 * (df + q + 1) * sum(log(diag(chol(A)))) 
  
  if (!isSymmetric(A)) {
    stop("A must be symmetric positive definite.")
  }
  
  exponent <- - 0.5 * sum(diag(S %*% chol2inv(chol(A))))
  
  if (log) {
    return(const + exponent)
  } else {
    return(exp(const + exponent))
  }
}

# log-posterior of the model parameters

log_posterior <- function(beta, Sigma, phi, nu, r, Y_vec, locations){
  
  prior_phi <- dunif(phi, 0, 1, log = TRUE)
  prior_nu <-  dlnorm(nu, meanlog = -1.2, sdlog = 1, log = TRUE)
  prior_Sigma <- dinvwish(A = Sigma,
                          df = q + 1,
                          S = diag(q),
                          log = TRUE)
  prior_beta <- dmatnorm(X = beta,
                         M = matrix(0, p, q),
                         U = (1e4)*diag(p),
                         V = Sigma,
                         log = TRUE)
  prior_r <- dunif(r, 0, 1, log = TRUE)
  
  
  return(log_likelihood(beta, Sigma, phi, nu, r, Y_vec, locations) +
           prior_beta +
           prior_Sigma +
           prior_phi +
           prior_nu +
           prior_r)
  
}

# Metropolis-Hastings algorithm for component-wise sampling

metro_gibbs <- function(beta, Sigma, phi, nu, r, 
                        niters, locations, Y_vec, X_vec,
                             # priors
                             M_prior ,
                             V_prior ,
                             S_prior,
                             df_prior,
                             #tuning params for MH
                             tuning_params) {
  
  beta_chain <- replicate(niters, matrix(0, p, q), simplify = F)
  Sigma_chain <- replicate(niters, matrix(0, q, q), simplify = F)
  phi_chain <- numeric(length = niters)
  nu_chain <- numeric(length = niters)
  r_chain <- numeric(length = niters)
  
  accept <- numeric(3)
  
  # Run Metropolis-Hastings
  for (iter in 2:niters) {
    
    if(iter %% ((niters)/10) == 0) print(paste0(100*(iter/(niters)), "%"))
    
    current_phi <- phi
    current_nu <- nu
    current_r <- r
    
    distmat <- rdist(locations)
    
    K.tilde <- r * Matern(distmat, 
                                  range = phi, 
                                  smoothness = nu) + (1- r)*diag(N)
    
    # Gibbs Update
    
    Y <- Y_vec
    dim(Y) <- c(N, q)
    
    X <- X_vec
    dim(X) <- c(N, p)
    
    K.tilde <- r * Matern(distmat, 
                          range = phi, 
                          smoothness = nu) + (1- r)*diag(N)
    
    prec.K.tilde <- chol2inv(chol(K.tilde))
    prec.V.prior <- chol2inv(chol(V_prior))

    # post_var of B
    V_tilde <-  chol2inv(chol(t(X) %*% prec.K.tilde %*% X + prec.V.prior))

    # post_mean of B
    M_tilde <- V_tilde %*% (t(X) %*% prec.K.tilde %*% Y + 
                              prec.V.prior %*% M_prior )
    
    
    
    resid <- (Y - X %*% beta)

    # post_scale matrix of Sigma
    S_tilde <- S_prior +  t(resid) %*% prec.K.tilde %*% resid 

    # post_df of Sigma
    df_tilde <- df_prior + N
    
    beta <- as.vector(vec(t(M_tilde))) +
      c(t(chol(Sigma %x% V_tilde))  %*% rnorm(p*q))
    
    dim(beta) = c(p,q)
    
    Sigma <- riwish(v = df_tilde, S = S_tilde)
    
    Sigma <- current_Sigma
    
    
    # Alternative: simulate n draws using rmniw function form mniw package
    
    # samples <- rMNIW(1, Lambda = M_tilde, Sigma = V_tilde, 
    #                  Psi = S_tilde, nu = df_tilde)
    # 
    # beta <- samples$X
    # Sigma <- samples$V
    
    
    
    # Metropolis update of the parameters
    
    current_phi <- rnorm(1, phi, sqrt(tuning_params[1]))
    
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
      
      accept[1] <- accept[1] + 1
      
    }
    
    current_nu <- rnorm(1, nu, sqrt(tuning_params[2]))
    
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
      
      accept[2] <- accept[2] + 1
      
    }
    
    current_r <- rnorm(1, r, sqrt(tuning_params[3]))
    
    
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
      accept[3] <- accept[3] + 1
      
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
beta <- true_beta
  
#matrix(rep(0, p*q), nrow = p, ncol = q)
Sigma <- true_Sigma
#diag(q)
phi <- 0.1
nu <- 0.5
r <- 0.8

# Number of iterations
niters <- 1e4

# Tuning parameters list

tuning_params <- c(5e-2, 8e-2, 8e-6)

# Run Metropolis-Hastings algorithm
theta_chain <- metro_gibbs(beta, Sigma, phi, nu, r, 
                                   niters, locations, Y_vec, X_vec,
                                   # priors
                                   M_prior = matrix(0, p, q) ,
                                   V_prior = 1e4*diag(p) ,
                                   S_prior = diag(q),
                                   df_prior = q + 1,
                                   #tuning params for MH
                                   tuning_params = tuning_params)

# Saving MCMC chain
list.save(theta_chain, file = "separable_mgp_metro_gibbs.Rdata")

theta_metro_gibbs <- list.load(file = "separable_mgp_metro_gibbs.Rdata")

# Traceplots

par(mfrow = c(1,1))
trace_phi <- plot.ts(theta_chain$phi_sample, ylab = "phi", main = "Traceplot of phi")
abline(h = true_phi, col = 'blue', lwd = 2)

trace_nu <- plot.ts(theta_chain$nu_sample, ylab = "nu", main = "Traceplot of nu")
abline(h = true_nu, col = 'blue', lwd = 2)

trace_r <- plot.ts(theta_chain$r_sample, ylab = "r", main = "Traceplot of r")
abline(h = true_r, col = 'blue', lwd = 2)

# acfplots

acf_phi <- acf(theta_chain$phi_sample, main = "ACF plot of phi", lag.max = 1e2)
acf_nu <- acf(theta_chain$nu_sample, main = "ACF plot of nu", lag.max = 1e2)
acf_r <- acf(theta_chain$r_sample, main = "ACF plot of r", lag.max = 1e2)

labels_beta <- list()
for (i in 1:p) {
  for (j in 1:q) {
    labels_beta <- c(labels_beta, 
                     list(paste0('beta (', i, ',', j, ')')))
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
        main = labels_beta[(i-1)*q + j], lag.max = 1e2)
    
  }
}


# Output Analysis of Sigma

labels_Sigma <- list()
for (i in 1:p) {
  for (j in 1:q) {
    labels_Sigma <- c(labels_Sigma, 
                      list(paste0('Sigma (', i, ',', j, ')')))
  }
}

par(mfrow=c(q,q))

for (i in 1:q) {
  for (j in 1:q) {
    plot(1:niters, sapply(theta_chain$Sigma_sample, function(x) x[i, j]), 
         type='l', col=1, xlab='Iteration', ylab='Sigma',
         main=labels_Sigma[(i-1)*q + j])
    abline(h = true_Sigma[i,j], col = 'blue', lwd = 2)
    
  }
}

par(mfrow=c(q,q))

for (i in 1:q) {
  for (j in 1:q) {
    
    acf(sapply(theta_chain$Sigma_sample, function(x) x[i, j]), 
        main=labels_Sigma[(i-1)*q + j],
        lag.max = 1e2)
    
  }
}
