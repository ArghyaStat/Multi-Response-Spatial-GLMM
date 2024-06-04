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
  
  
  denom <- -0.5 * (p * q * log(2 * pi) + q * sum(log(diag(chol.V))) + p * sum(log(diag(chol.U))))
  
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
  
  const <- (0.5 * q * df) * sum(log(diag(chol(S)))) 
  - (0.5 * df) * log(2) - lgamma(0.5 * q * df )
  
  if (!isSymmetric(A)) {
    stop("A must be symmetric positive definite.")
  }
  
  exponent <- - 0.5 * (df + q + 1) * sum(log(diag(chol(A)))) 
  - 0.5 * sum(diag(S %*% chol2inv(chol(A))))
  
  if (log) {
    return(const + exponent)
  } else {
    return(exp(const + exponent))
  }
}


# MCMC Set up
# Conditional Log-likelihood function for W.tilde | Y using Cholesky decomposition

log.likelihood <- function(W.tilde, Y) {
    
    lambda1 <- exp(W.tilde[,1])
    lambda2 <- exp(W.tilde[,2])
    
    
    like.count1 <- sum(dpois(Y[,1], lambda = lambda1, log = T))
    like.count2 <- sum(dpois(Y[,2], lambda = lambda2, log = T))  
    
    
    like.out <- like.count1 + like.count2
  
  return(like.out)
}


# log-posterior of the model parameters

log.posterior <- function(W.tilde, beta, Sigma, phi, nu, r, Y, X, locations){
  
  
  prior.Sigma <- dinvwish(A = Sigma, 
                          df = q + 1, 
                          S = diag(q), 
                          log = TRUE)
  prior.beta <- dmatnorm(X = beta, 
                         M = matrix(0, p, q), 
                         U = 1e4 * diag(p),
                         V = Sigma,
                         log = TRUE)
  
  distmat <- rdist(locations)
  
  if(phi > 0 & nu > 0 & r > 0 & r < 1){
    
  K.tilde <- r*Matern(distmat, range = phi, smoothness = nu) + (1- r)*diag(N)
  prior.W.tilde <- dmatnorm(X = W.tilde, 
                            M = X %*% beta, 
                            U = K.tilde,
                            V = Sigma,
                            log = TRUE)
  
  prior.phi <- dunif(phi, 0, 1, log = TRUE)
  prior.nu <-  dlnorm(nu, meanlog = -1.2, sdlog = 1, log = TRUE)
  prior.r <- dunif(r, 0, 1, log = TRUE)
  
  return(log.likelihood(W.tilde, Y) +
           prior.W.tilde + 
           prior.beta +
           prior.Sigma +
           prior.phi +
           prior.nu +
           prior.r)
  }else{
    
    return(-Inf)
  }
  
  
}
