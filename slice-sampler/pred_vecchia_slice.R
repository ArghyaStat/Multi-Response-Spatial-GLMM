
W.pred.ord <- function(W.obs.ord, beta, chol.Sigma, phi, nu, m,
                       X.obs.ord, X.pred.ord, U.joint, U.pred.col, U.obs.col, 
                       chol.W.pred.var, N.obs, N.pred, q, p){

  
  W.pred.mean.part <- U.obs.col %*% (W.obs.ord - X.obs.ord %*% beta)
  W.pred.mean.tilde <- crossprod.spam(U.pred.col,  W.pred.mean.part)
 
                                      
  W.pred.samp <- (matrix(rnorm(N.pred*q, 0, 1), nrow = N.pred , ncol = q) %*% chol.Sigma)
  
  W.pred.temp <- chol.W.pred.var %*% W.pred.mean.tilde +  W.pred.samp 
  
  W.pred.ord <- X.pred.ord %*% beta -  crossprod(chol.W.pred.var,  W.pred.temp)
  
  W.pred.ord
  
  
}


Y.pred.ord <- function(W.pred.ord, N.pred, q, family){
  
  # Validate inputs
  if (length(family) != q) {
    stop("The length of 'family' must be equal to q.")
  }
  
  Y.mat <- matrix(NA, nrow = N.pred, ncol = q)
  
  # Validate inputs
  if (length(family) != q) {
    stop("The length of 'family' must be equal to q.")
  }
  
  
  # Generate responses based on family
  for (j in 1:q) {
    if (family[j] == "Gaussian") {
      Y.mat[, j] <- rnorm(N.pred, mean = W.pred.ord[, j], sd = 1)
    } else if (family[j] == "Poisson") {
      Y.mat[, j] <- rpois(N.pred, lambda = exp(W.pred.ord[, j]))
    } else if (family[j] == "Binomial") {
      Y.mat[, j] <- rbinom(N.pred, size = 1, prob = (exp(-W.pred.ord[, j]) / (1 + exp(-W.pred.ord[, j]))))
    } else if (family[j] == "Gamma") {
      Y.mat[, j] <- rgamma(N.pred, shape = exp(W.pred.ord[, j]), rate = 1)
    } else if (family[j] == "Negative-Binomial") {
      Y.mat[, j] <- MASS::rnegbin(N.pred, mu = exp(W.pred.ord[, j]), theta = 1) # Theta as dispersion parameter
    } else {
      stop(paste("Unsupported family:", family[j]))
    }
  }
  
  Y.mat
}


# thinning <- function(pred.iters, niters, burnin,
#                      W.ord.samples, beta.samples, 
#                      Sigma.samples, chol.Sigma.samples, phi.samples){
#   
#   
#   
#   if ((niters-burnin) %% pred.iters != 0) {
#     stop("niters must be divisible by pred.iters")
#   }
#   
#   thin.block <- (niters-burnin)/pred.iters
#   thinned.indices <- seq(burnin + 1, niters, by = thin.block)
#   
#   thinned.W.ord.samples <- list()
#   thinned.beta.samples <- list()
#   thinned.Sigma.samples <- list()
#   thinned.chol.Sigma.samples <- list()
#   
#   
#   for (i in thinned.indices) {
#     
#     thinned.W.ord.samples <- c(thinned.W.ord.samples, list(W.ord.samples[[i]]))
#     thinned.beta.samples <- c(thinned.beta.samples, list(beta.samples[[i]]))
#     thinned.Sigma.samples <- c(thinned.Sigma.samples, list(Sigma.samples[[i]]))
#     thinned.chol.Sigma.samples <- c(thinned.chol.Sigma.samples, list(chol.Sigma.samples[[i]]))
#     
#   }
#   
#   thinned.phi.samples <- phi.samples[thinned.indices]
#   
#   thinned.samples <- list("thinned.W.ord.samples" = thinned.W.ord.samples,
#                           "thinned.beta.samples" = thinned.beta.samples,
#                           "thinned.Sigma.samples" = thinned.Sigma.samples,
#                           "thinned.chol.Sigma.samples" = thinned.chol.Sigma.samples,
#                           "thinned.phi.samples" = thinned.phi.samples
#   )
#   
#   return(thinned.samples)
#   
# }


predictive.samples <- function(W.obs.ord.samples, beta.samples, chol.Sigma.samples, 
                               phi.samples, nu, m,
                               X.obs.ord, X.pred.ord, U.joint, U.pred.col, U.obs.col, 
                               chol.W.pred.var, N.obs, N.pred, pred.iters, q, p){
  
  Y.pred.ord <- replicate(pred.iters, matrix(0, N.pred, q), simplify = F)
  W.pred.ord <- replicate(pred.iters, matrix(0, N.pred, q), simplify = F)
  
  for(i in 1:pred.iters){
    
    W.pred.ord[[i]] <-  W.pred.ord(W.obs.ord = W.obs.ord.samples[[i]], 
                                   beta =  beta.samples[[i]],
                                   chol.Sigma = chol.Sigma.samples[[i]], 
                                   phi = phi.samples[[i]], nu, m,
                                   X.obs.ord, X.pred.ord, U.joint, U.pred.col, U.obs.col, 
                                   chol.W.pred.var, N.obs, N.pred, q, p)
    Y.pred.ord[[i]] <- Y.pred.ord(W.pred.ord[[i]], N.pred, q, family)
  }
  
  Y.pred.ord
  
}