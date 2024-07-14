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
  
  
  denom <- -0.5 * (p * q * log(2 * pi) + 2*q * sum(log(diag(chol.V))) + 2*p * sum(log(diag(chol.U))))
  
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

log.likelihood <- function(W.tilde, Y, N) {
    
    lambda1 <- exp(W.tilde[,1])
    lambda2 <- exp(W.tilde[,2])
    
    
    like.count1 <- sum(dpois(Y[,1], lambda = lambda1, log = T))
    like.count2 <- sum(dpois(Y[,2], lambda = lambda2, log = T))  
    
    
    like.out <- like.count1 + like.count2
  
  return(like.out)
}



cormat.update <- function(distmat, phi, nu, r, N){
  
  # required library(fields)
  
  K.tilde <- r* Matern(distmat, range = phi, smoothness = nu) + (1-r)* diag(N)
  
  chol.K.tilde <- chol(K.tilde)
  K.tilde.inv <- chol2inv(chol.K.tilde)
  
  list(cormat = K.tilde,
       cormat.inv = K.tilde.inv)
}


update.W.tilde <- function(W.tilde, beta, Sigma, r, phi, nu, Y, X, N, 
                           K.tilde, distmat, acc.W.tilde, tuning.W.tilde){
  
  q <- nrow(Sigma)
  
  post.W.tilde <- log.likelihood(W.tilde, Y, N) +
                            dmatnorm(W.tilde,
                                     M = X %*% beta,
                                     U = K.tilde,
                                     V = Sigma,
                                     log = TRUE)
  W.error <- sqrt(tuning.W.tilde) * matrix(rnorm(N*q, 0, 1), nrow = N, ncol = q)
  
  # Sampling form matrix-normal using cholesky decomposition
  
  can.W.tilde <- W.tilde + W.error
  
  post.can.W.tilde <- log.likelihood(can.W.tilde, Y, N) +
                            dmatnorm(can.W.tilde,
                                     M = X %*% beta,
                                     U = K.tilde,
                                     V = Sigma,
                                    log = TRUE)
  
  log.r.W.tilde <- post.can.W.tilde - post.W.tilde
  
  if(log(runif(1)) < log.r.W.tilde){
    
    W.tilde <- can.W.tilde
    acc.W.tilde <- acc.W.tilde + 1
    
  }
  
  results = list(W.tilde = W.tilde,  acc.W.tilde = acc.W.tilde)
  
  results
}



update.beta <- function(W.tilde, Sigma, X, K.tilde.inv, M.prior, V.prior){
  
  chol.V.prior <- chol(V.prior)
  prec.V.prior <- chol2inv(chol.V.prior)
  
  chol.Sigma <- chol(Sigma)
  
  V.tilde <- chol2inv(chol(t(X) %*% K.tilde.inv %*% X + prec.V.prior))
  M.tilde <- V.tilde %*% (t(X) %*% K.tilde.inv %*% W.tilde + prec.V.prior %*% M.prior)
  
  chol.V.tilde <- chol(V.tilde)
  
  beta.samp <- matrix(rnorm(p*q, 0, 1), nrow = p , ncol = q)
  
  beta <- M.tilde + t(chol.V.tilde) %*% beta.samp %*% chol.Sigma
  
  beta
  
}
  

update.Sigma <- function(W.tilde, beta, X, K.tilde.inv, 
                         M.prior, V.prior, S.prior, df.prior, N){
  
  p <- nrow(V.prior)
  N <- nrow(W.tilde)
  
  chol.V.prior <- chol(V.prior)
  prec.V.prior <- chol2inv(chol.V.prior)
  
  resid.W.tilde <- W.tilde - X %*% beta
  rss.W.tilde <- t(resid.W.tilde) %*% K.tilde.inv %*% resid.W.tilde
  
  rss.beta <- t(beta - M.prior) %*% prec.V.prior %*% (beta - M.prior)
  
  S.tilde <- S.prior + rss.W.tilde + rss.beta
  df.tilde <- df.prior + N + p
    
  Sigma <- riwish(v = df.tilde, S = S.tilde)
  
  Sigma
}


update.phi <- function(phi, beta, Sigma, nu, r, W.tilde, distmat, 
                       K.tilde, K.tilde.inv, acc.phi, tuning.phi){
  
  post.phi <- dmatnorm(W.tilde,
                       M = X %*% beta,
                       U = K.tilde,
                       V = Sigma,
                       log = TRUE) + dunif(phi, 0, 1, log = TRUE)
  
  
  can.phi <- rnorm(1, phi, sqrt(tuning.phi))
  
  # cormat details for candidate phi
  
  if(can.phi > 0 ){
  
  can.cormat.details <- cormat.update(distmat, can.phi, nu, r, N)
  
  can.K.tilde <- can.cormat.details$cormat
  can.K.tilde.inv <- can.cormat.details$cormat.inv
  
  
  post.can.phi <- dmatnorm(W.tilde,
                           M = X %*% beta,
                           U = can.K.tilde,
                           V = Sigma,
                           log = TRUE) + dunif(can.phi, 0, 1, log = TRUE)
  
  
  
  log.r.phi <- post.can.phi - post.phi
  
  if(log(runif(1)) < log.r.phi){
    
    phi = can.phi
    K.tilde = can.K.tilde
    K.tilde.inv = can.K.tilde.inv
    
    acc.phi = acc.phi + 1
  }
 }
  
  results = list(phi = phi, cormat = K.tilde, 
                 cormat.inv = K.tilde.inv, acc.phi = acc.phi)
  
  results
  
}



update.nu <- function(phi, beta, Sigma, nu, r, W.tilde, distmat, 
                      K.tilde, K.tilde.inv, acc.nu, tuning.nu){
  
  post.nu <- dmatnorm(W.tilde,
                       M = X %*% beta,
                       U = K.tilde,
                       V = Sigma,
                       log = TRUE) + 
             dlnorm(nu, meanlog = - 1.2, sdlog = 1, log = TRUE)
  
  
    can.nu <- rnorm(1, nu, sqrt(tuning.nu))
  
  # cormat details for candidate phi
  
  if(can.nu > 0){
    
    can.cormat.details <- cormat.update(distmat, phi, can.nu, r, N)
    
    can.K.tilde <- can.cormat.details$cormat
    can.K.tilde.inv <- can.cormat.details$cormat.inv
    
    
    post.can.nu <- dmatnorm(W.tilde,
                             M = X %*% beta,
                             U = can.K.tilde,
                             V = Sigma,
                             log = TRUE) + 
                   dlnorm(can.nu, meanlog = -1.2, sdlog = 1, log = TRUE)
    
    
    
    log.r.nu <- post.can.nu - post.nu
    
    if(log(runif(1)) < log.r.nu){
      
      nu = can.nu
      K.tilde = can.K.tilde
      K.tilde.inv = can.K.tilde.inv
      
      acc.nu = acc.nu + 1
    }
  }
  
  results = list(nu = nu, cormat = K.tilde, 
                 cormat.inv = K.tilde.inv, acc.nu = acc.nu)
  
  results
  
}



update.r <- function(phi, beta, Sigma, nu, r, W.tilde, distmat, 
                     K.tilde, K.tilde.inv, acc.r, tuning.r){
  
  post.r <- dmatnorm(W.tilde,
                      M = X %*% beta,
                      U = K.tilde,
                      V = Sigma,
                      log = TRUE) + 
            dunif(r, 0, 1, log = TRUE)
  
  
  can.r <- rnorm(1, r, sqrt(tuning.r))
  
  # cormat details for candidate phi
  
  if(can.r > 0 & can.r < 1){
    
    can.cormat.details <- cormat.update(distmat, phi, nu, can.r, N)
    
    can.K.tilde <- can.cormat.details$cormat
    can.K.tilde.inv <- can.cormat.details$cormat.inv
    
    
    post.can.r <- dmatnorm(W.tilde,
                             M = X %*% beta,
                             U = can.K.tilde,
                             V = Sigma,
                             log = TRUE) + 
                    dunif(can.r, 0, 1, log = TRUE)
    
    
    
    log.r.r <- post.can.r - post.r
    
    if(log(runif(1)) < log.r.r){
      
      r = can.r
      K.tilde = can.K.tilde
      K.tilde.inv = can.K.tilde.inv
      
      acc.r = acc.r + 1
    }
  }
  
  results = list(r = r, cormat = K.tilde, 
                 cormat.inv = K.tilde.inv, acc.r = acc.r)
  
  results
  
}

