# Define own function for matrix-normal density
dmatnorm <- function(X, M, U, V, log = FALSE) {

  
  if (!isSymmetric(V)) {
    stop("V must be symmetric positive definite.")
  }
  if (!isSymmetric(U)) {
    stop("U must be symmetric positive definite.")
  }
  
  p <- nrow(X)
  q <- ncol(X)
  
  chol.U <- chol(U)
  chol.V <- chol(V)
  
  
  denom <- -0.5 * (p * q * log(2 * pi) + 2*p * sum(log(diag(chol.V))) + 2*q * sum(log(diag(chol.U))))
  
  V.prec <- chol2inv(chol.V)
  U.prec <- chol2inv(chol.U)
  
  exponent <- -0.5 * sum(diag(V.prec %*% t(X - M) %*% U.prec %*% (X-M)))
  
  if (log) {
    return(denom + exponent)
  } else {
    return(exp(denom + exponent))
  }
}


# Define function for multivariate-normal density

dmvn <- function(x, mu.x, cov.x, log = FALSE) {
  
  
  if (!isSymmetric(cov.x)) {
    stop("cov.x must be symmetric positive definite.")
  }
  
  k <- length(x)
  chol.cov.x <- chol(cov.x)
  prec.x <- chol2inv(chol.cov.x)
  
  denom <- -0.5 * (k * log(2 * pi) + 2* sum(log(diag(chol.cov.x))))
  
  exponent <- -0.5 * t(x - mu.x) %*% prec.x %*% (x - mu.x)
  
  if (log) {
    return(denom + exponent)
  } else {
    return(exp(denom + exponent))
  }
}


# MCMC Set up
# Conditional Log-likelihood function for W.tilde | Y using Cholesky decomposition

log.likelihood <- function(W.tilde.obs, Y.obs) {

    N.obs <- nrow(W.tilde.obs)

    like.count <-  0

    like.cont <- 0


    for(i in 1:N.obs){

      like.cont <-  like.cont + dnorm(Y.obs[i,1],  mean = W.tilde.obs[i,1],  1, log = T)
      
      like.count <- like.count + dpois(Y.obs[i,2], lambda = exp(W.tilde.obs[i,2]), log = T)

    }

    like.out <- like.count + like.cont

  return(like.out)
}



cormat.update <- function(distmat.obs, phi, nu, r, N.obs){
  
  # required library(fields)
  
  K.tilde.obs <- r* Matern(distmat.obs, range = phi, smoothness = nu) + (1-r)* diag(N.obs)
  
  chol.K.tilde.obs <- chol(K.tilde.obs)
  K.tilde.obs.inv <- chol2inv(chol.K.tilde.obs)
  
  list(cormat.obs = K.tilde.obs,
       cormat.obs.inv = K.tilde.obs.inv)
}


# Independent Gaussian proposal with conditional MAP as location

update.W.tilde <- function(W.tilde.obs, beta, Sigma, K.tilde.obs, K.tilde.obs.inv,
                           Y.obs, X.obs, N.obs, distmat.obs, H.inv,
                           W.tilde.mle, acc.W.tilde, tuning.W.tilde){
  
  # vectorize both Y and W from matrices
  W.tilde.obs.vec <- vec(W.tilde.obs)
  Y.vec <- vec(Y.obs)
  
  # Prior covariance of W
  Omega <- kronecker(Sigma, K.tilde.obs)
  
  Sigma.inv <- chol2inv(chol(Sigma))
  Omega.inv <- kronecker(Sigma.inv, K.tilde.obs.inv)
  
  # Finding hessian of MAP estimator 
  Omega.map <- chol2inv(chol(Omega.inv + H.inv))
  
  # Finding the conditional MAP estimator
  
  W.tilde.map <- Omega.map %*% (Omega.inv %*% vec(X.obs %*% beta) + H.inv %*% W.tilde.mle)
  
  chol.Omega.map <- chol(Omega.map)
  
  # Calculation of target density at correct W
  
  post.W.tilde <- log.likelihood(W.tilde.obs, Y.obs) +
                  dmvn(x = c(W.tilde.obs.vec),
                  mu.x = c(vec(X.obs %*% beta)),
                  cov.x = Omega, log = TRUE)
  
  # Generating candidate W
  
  W.error <- sqrt(tuning.W.tilde) * c(t(chol.Omega.map) %*% rnorm(N.obs*q))
  can.W.tilde.obs.vec <- W.tilde.obs.vec + W.error
  
  
  # calculating target at candidate W 
  
  can.W.tilde.obs <- matrix(can.W.tilde.obs.vec, nrow = N.obs, ncol = q)
  
  post.can.W.tilde <- log.likelihood(can.W.tilde.obs, Y.obs) +
                      dmvn(x = c(can.W.tilde.obs.vec),
                      mu.x = c(vec(X.obs %*% beta)),
                      cov.x = Omega, log = TRUE)
  
  # q(candidate | current)
  fwd.density.W.tilde <- dmvn(x = c(can.W.tilde.obs.vec),
                              mu.x = c(W.tilde.map),
                              cov.x = sqrt(tuning.W.tilde) * Omega.map, log = TRUE)
  # q(current | candidate)
  rev.density.W.tilde <- dmvn(x = c(W.tilde.obs.vec),
                              mu.x = c(W.tilde.map),
                              cov.x = sqrt(tuning.W.tilde) * Omega.map, log = TRUE)
  
  # Calculating MH acceptance ratio in log-scale  
  log.r.W.tilde <- post.can.W.tilde - post.W.tilde + rev.density.W.tilde - fwd.density.W.tilde
  
  if(log(runif(1)) < log.r.W.tilde){
    
    W.tilde.obs.vec <- can.W.tilde.obs.vec
    acc.W.tilde <- acc.W.tilde + 1
    
  }
  
  # reconverting the candidate from vec to matrix
  
  W.tilde.obs = matrix(W.tilde.obs.vec, nrow = N.obs, ncol = q)
  
  results = list(W.tilde.obs = W.tilde.obs,  acc.W.tilde = acc.W.tilde)
  
  results
}


# update of phi with truncated normal proposal centered at current phi

update.phi <- function(phi, beta, Sigma, nu, r, W.tilde.obs, distmat.obs, 
                       K.tilde.obs, K.tilde.obs.inv, acc.phi, tuning.phi, b_phi){
  
  # Calculate the posterior for the current phi
  post.phi <- dmatnorm(W.tilde.obs,
                       M = X.obs %*% beta,
                       U = K.tilde.obs,
                       V = Sigma,
                       log = TRUE) + dunif(phi, 0, b_phi, log = TRUE)
  
  # Generate a candidate phi using a truncated Gaussian distribution
  can.phi <- rtruncnorm(1, a = 0, b = b_phi, mean = phi, sd = sqrt(tuning.phi))
  
  if(can.phi > 0 && can.phi < b_phi){
    
    # Update the correlation matrix details for the candidate phi
    can.cormat.details <- cormat.update(distmat.obs, can.phi, nu, r, N.obs)
    can.K.tilde.obs <- can.cormat.details$cormat.obs
    can.K.tilde.obs.inv <- can.cormat.details$cormat.obs.inv
    
    # Calculate the posterior for the candidate phi
    post.can.phi <- dmatnorm(W.tilde.obs,
                             M = X.obs %*% beta,
                             U = can.K.tilde.obs,
                             V = Sigma,
                             log = TRUE) + dunif(can.phi, 0, b_phi, log = TRUE)
    
    # Compute the truncated normal density at the current and candidate phi
    q_curr_given_can <- dtruncnorm(phi, a = 0, b = b_phi, mean = can.phi, sd = sqrt(tuning.phi))
    q_can_given_curr <- dtruncnorm(can.phi, a = 0, b = b_phi, mean = phi, sd = sqrt(tuning.phi))
    
    # Calculate the log acceptance ratio with correction for asymmetry
    log.r.phi <- (post.can.phi - post.phi) + 
      log(q_curr_given_can) - log(q_can_given_curr)
    
    # Metropolis-Hastings acceptance step
    if(log(runif(1)) < log.r.phi){
      phi <- can.phi
      K.tilde.obs <- can.K.tilde.obs
      K.tilde.obs.inv <- can.K.tilde.obs.inv
      acc.phi <- acc.phi + 1
    }
  }
  
  results <- list(phi = phi, cormat.obs = K.tilde.obs, 
                  cormat.obs.inv = K.tilde.obs.inv, acc.phi = acc.phi)
  
  return(results)
}



update.beta <- function(W.tilde.obs, Sigma, X.obs, K.tilde.obs.inv, M.tilde, chol.V.tilde){
  
  
  chol.Sigma <- chol(Sigma)
  
  beta.samp <- matrix(rnorm(p*q, 0, 1), nrow = p , ncol = q)
  
  beta <- M.tilde + t(chol.V.tilde) %*% beta.samp %*% chol.Sigma
  
  beta

  
}
  

update.Sigma <- function(W.tilde.obs, X.obs, K.tilde.obs.inv, 
                         M.prior, prec.V.prior, M.tilde, prec.V.tilde, 
                         S.prior, df.prior, N.obs){
  
  
  rss.W.tilde.obs <- t(W.tilde.obs) %*% K.tilde.obs.inv %*% W.tilde.obs
  
  rss.beta.prior <- t(M.prior) %*% prec.V.prior %*% (M.prior)
  rss.beta.posterior <- t(M.tilde) %*% prec.V.tilde %*% (M.tilde)
  
  S.tilde <- S.prior + rss.W.tilde.obs + rss.beta.prior - rss.beta.posterior
  df.tilde <- df.prior + N.obs 
    
  Sigma <- riwish(v = df.tilde, S = S.tilde)
  
  Sigma
}





update.nu <- function(phi, beta, Sigma, nu, r, W.tilde.obs, distmat.obs, 
                      K.tilde.obs, K.tilde.obs.inv, acc.nu, tuning.nu){
  
  post.nu <- dmatnorm(W.tilde.obs,
                       M = X.obs %*% beta,
                       U = K.tilde.obs,
                       V = Sigma,
                       log = TRUE) + 
             dlnorm(nu, meanlog = -log(2) - 0.5, sdlog = 1, log = TRUE)
  
  
    can.nu <- rnorm(1, nu, sqrt(tuning.nu))
  
  # cormat details for candidate nu
  
  if(can.nu > 0){
    
    can.cormat.details <- cormat.update(distmat.obs, phi, can.nu, r, N.obs)
    
    can.K.tilde.obs <- can.cormat.details$cormat.obs
    can.K.tilde.obs.inv <- can.cormat.details$cormat.obs.inv
    
    
    post.can.nu <- dmatnorm(W.tilde.obs,
                             M = X.obs %*% beta,
                             U = can.K.tilde.obs,
                             V = Sigma,
                             log = TRUE) + 
                   dlnorm(can.nu,  meanlog = -log(2) - 0.5, sdlog = 1, log = TRUE)
    
    
    
    log.r.nu <- post.can.nu - post.nu
    
    if(log(runif(1)) < log.r.nu){
      
      nu = can.nu
      K.tilde.obs = can.K.tilde.obs
      K.tilde.obs.inv = can.K.tilde.obs.inv
      
      acc.nu = acc.nu + 1
    }
  }
  
  results = list(nu = nu, cormat.obs = K.tilde.obs, 
                 cormat.obs.inv = K.tilde.obs.inv, acc.nu = acc.nu)
  
  results
  
}



update.r <- function(phi, beta, Sigma, nu, r, W.tilde.obs, distmat, 
                     K.tilde.obs, K.tilde.obs.inv, acc.r, tuning.r){
  
  post.r <- dmatnorm(W.tilde.obs,
                      M = X.obs %*% beta,
                      U = K.tilde.obs,
                      V = Sigma,
                      log = TRUE) + 
            dunif(r, 0, 1, log = TRUE)
  
  
  can.r <- rnorm(1, r, sqrt(tuning.r))
  
  # cormat details for candidate phi
  
  if(can.r > 0 & can.r < 1){
    
    can.cormat.details <- cormat.update(distmat.obs, phi, nu, can.r, N.obs)
    
    can.K.tilde.obs <- can.cormat.details$cormat.obs
    can.K.tilde.obs.inv <- can.cormat.details$cormat.obs.inv
    
    
    post.can.r <- dmatnorm(W.tilde.obs,
                             M = X.obs %*% beta,
                             U = can.K.tilde.obs,
                             V = Sigma,
                             log = TRUE) + 
                    dunif(can.r, 0, 1, log = TRUE)
    
    
    
    log.r.r <- post.can.r - post.r
    
    if(log(runif(1)) < log.r.r){
      
      r = can.r
      K.tilde.obs = can.K.tilde.obs
      K.tilde.obs.inv = can.K.tilde.obs.inv
      
      acc.r = acc.r + 1
    }
  }
  
  results = list(r = r, cormat.obs = K.tilde.obs, 
                 cormat.obs.inv = K.tilde.obs.inv, acc.r = acc.r)
  
  results
  
}



# Proposal 1

# update.W.tilde <- function(W.tilde.obs, beta, Sigma, K.tilde.obs, K.tilde.obs.inv, 
#                            Y.obs, X.obs, N.obs, distmat.obs, H.inv,
#                            W.tilde.mle, acc.W.tilde, tuning.W.tilde){
#   
#   W.tilde.obs.vec <- vec(W.tilde.obs)
#   Y.vec <- vec(Y.obs)
#   
#   Omega <- kronecker(Sigma, K.tilde.obs)
#   
#   Sigma.inv <- chol2inv(chol(Sigma))
#   Omega.inv <- kronecker(Sigma.inv, K.tilde.obs.inv)
#   
#   Omega.map <- chol2inv(chol(Omega.inv + H.inv))
#   W.tilde.map <- Omega.map %*% (Omega.inv %*% vec(X.obs %*% beta) + H.inv %*% W.tilde.mle)
#   
#   chol.Omega.map <- chol(Omega.map)
#   
#   post.W.tilde <- log.likelihood(W.tilde.obs, Y.obs) + 
#                   dmatnorm(W.tilde.obs,
#                            M = X.obs %*% beta,
#                            U = K.tilde.obs,
#                            V = Sigma,
#                            log = TRUE)
#   
#   W.tilde.mean.vec <- W.tilde.map + sqrt(1 - tuning.W.tilde)*(W.tilde.obs.vec - W.tilde.map)
#   
#   
#   W.error <- sqrt(tuning.W.tilde) * c(t(chol.Omega.map) %*% rnorm(N.obs*q))  
#   
#   can.W.tilde.obs.vec <- W.tilde.mean.vec + W.error
#   
#   can.W.tilde.mean.vec <- W.tilde.map + sqrt(1 - tuning.W.tilde)*(can.W.tilde.obs.vec - W.tilde.map)
#   
#   can.W.tilde.obs <- matrix(can.W.tilde.obs.vec, nrow = N.obs, ncol = q)
#   
#   post.can.W.tilde <- log.likelihood(can.W.tilde.obs, Y.obs) +
#     dmatnorm(can.W.tilde.obs,
#              M = X.obs %*% beta,
#              U = K.tilde.obs,
#              V = Sigma,
#              log = TRUE)
# 
#   fwd.density.W.tilde <- dmvnorm(x = can.W.tilde.obs.vec,
#                                  mu.x = W.tilde.mean.vec, 
#                                  cov.x = sqrt(tuning.W.tilde)*Omega.map, 
#                                  log = TRUE)
# 
#   rev.density.W.tilde <- dmvnorm(x = W.tilde.obs.vec,
#                          mu.x = can.W.tilde.mean.vec, 
#                          cov.x = sqrt(tuning.W.tilde)*Omega.map, 
#                          log = TRUE)
# 
#   
#   
#   
#   log.r.W.tilde <- post.can.W.tilde - post.W.tilde + 
#                    rev.density.W.tilde - fwd.density.W.tilde
#   
#   if(log(runif(1)) < log.r.W.tilde){
#     
#     W.tilde.obs <- can.W.tilde.obs
#     acc.W.tilde <- acc.W.tilde + 1
#     
#   }
#   
#   results = list(W.tilde.obs = W.tilde.obs,  acc.W.tilde = acc.W.tilde)
#   
#   results
# }

# update.phi <- function(phi, beta, Sigma, nu, r, W.tilde.obs, distmat.obs, 
#                        K.tilde.obs, K.tilde.obs.inv, acc.phi, tuning.phi){
#   
#   post.phi <- dmatnorm(W.tilde.obs,
#                        M = X.obs %*% beta,
#                        U = K.tilde.obs,
#                        V = Sigma,
#                        log = TRUE) + dunif(phi, 0, diameter/log(20), log = TRUE)
#     
#     
#   can.phi <- rnorm(1, phi, sqrt(tuning.phi))
#   
#   # cormat details for candidate phi
#   
#   if(can.phi > 0){
#   
#   can.cormat.details <- cormat.update(distmat.obs, can.phi, nu, r, N.obs)
#   
#   can.K.tilde.obs <- can.cormat.details$cormat.obs
#   can.K.tilde.obs.inv <- can.cormat.details$cormat.obs.inv
#   
#   
#   post.can.phi <- dmatnorm(W.tilde.obs,
#                            M = X.obs %*% beta,
#                            U = can.K.tilde.obs,
#                            V = Sigma,
#                            log = TRUE) + dunif(can.phi, 0, diameter/log(20), log = TRUE)
#   
#   
#   
#   log.r.phi <- post.can.phi - post.phi
#   
#   if(log(runif(1)) < log.r.phi){
#     
#     phi = can.phi
#     K.tilde.obs = can.K.tilde.obs
#     K.tilde.obs.inv = can.K.tilde.obs.inv
#     
#     acc.phi = acc.phi + 1
#   }
#  }
#   
#   results = list(phi = phi, cormat.obs = K.tilde.obs, 
#                  cormat.obs.inv = K.tilde.obs.inv, acc.phi = acc.phi)
#   
#   results
#   
# }


