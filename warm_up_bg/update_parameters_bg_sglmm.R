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


# Define own function for inverse-Wishart density

dinvwish <- function(A, df, S, log = FALSE) {
  
  q <- nrow(A)
  
  const <- (q * df) * sum(log(diag(chol(S)))) 
  - (0.5 * df) * log(2) - lgamma(0.5 * q * df )
  
  if (!isSymmetric(A)) {
    stop("A must be symmetric positive definite.")
  }
  
  exponent <- - (df + q + 1) * sum(log(diag(chol(A)))) 
  - 0.5 * sum(diag(S %*% chol2inv(chol(A))))
  
  if (log) {
    return(const + exponent)
  } else {
    return(exp(const + exponent))
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


update.W.tilde <- function(W.tilde.obs, beta, Sigma, K.tilde.obs,
                           Y.obs, X.obs, N.obs, distmat.obs, 
                           acc.W.tilde, var.W.tilde, tuning.W.tilde){
  
  q <- nrow(Sigma)
  
  post.W.tilde <- log.likelihood(W.tilde.obs, Y.obs) +
                            dmatnorm(W.tilde.obs,
                                     M = X.obs %*% beta,
                                     U = K.tilde.obs,
                                     V = Sigma,
                                     log = TRUE)
  W.error <- sqrt(tuning.W.tilde) * matrix(rnorm(N.obs * q, mean = 0, sd = sqrt(var.W.tilde)), 
                              nrow = N.obs, ncol = q)

  # Sampling form matrix-normal using cholesky decomposition
  
  can.W.tilde.obs <- W.tilde.obs + W.error
  
  post.can.W.tilde <- log.likelihood(can.W.tilde.obs, Y.obs) +
                            dmatnorm(can.W.tilde.obs,
                                     M = X.obs %*% beta,
                                     U = K.tilde.obs,
                                     V = Sigma,
                                     log = TRUE)
  
  log.r.W.tilde <- post.can.W.tilde - post.W.tilde
  
  if(log(runif(1)) < log.r.W.tilde){
    
    W.tilde.obs <- can.W.tilde.obs
    acc.W.tilde <- acc.W.tilde + 1
    
  }
  
  results = list(W.tilde.obs = W.tilde.obs,  acc.W.tilde = acc.W.tilde)
  
  results
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


update.phi <- function(phi, beta, Sigma, nu, r, W.tilde.obs, distmat.obs, 
                       K.tilde.obs, K.tilde.obs.inv, acc.phi, tuning.phi){
  
  post.phi <- dmatnorm(W.tilde.obs,
                       M = X.obs %*% beta,
                       U = K.tilde.obs,
                       V = Sigma,
                       log = TRUE) + dunif(phi, 0, diameter/log(20), log = TRUE)
    
    
  can.phi <- rnorm(1, phi, sqrt(tuning.phi))
  
  # cormat details for candidate phi
  
  if(can.phi > 0){
  
  can.cormat.details <- cormat.update(distmat.obs, can.phi, nu, r, N.obs)
  
  can.K.tilde.obs <- can.cormat.details$cormat.obs
  can.K.tilde.obs.inv <- can.cormat.details$cormat.obs.inv
  
  
  post.can.phi <- dmatnorm(W.tilde.obs,
                           M = X.obs %*% beta,
                           U = can.K.tilde.obs,
                           V = Sigma,
                           log = TRUE) + dunif(can.phi, 0, diameter/log(20), log = TRUE)
  
  
  
  log.r.phi <- post.can.phi - post.phi
  
  if(log(runif(1)) < log.r.phi){
    
    phi = can.phi
    K.tilde.obs = can.K.tilde.obs
    K.tilde.obs.inv = can.K.tilde.obs.inv
    
    acc.phi = acc.phi + 1
  }
 }
  
  results = list(phi = phi, cormat.obs = K.tilde.obs, 
                 cormat.obs.inv = K.tilde.obs.inv, acc.phi = acc.phi)
  
  results
  
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

