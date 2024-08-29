# mcmc

mcmc.sglmm <- function(Y.obs, X.obs, obs.loc,
                       #parameters
                       W.tilde.obs, beta, Sigma, phi, nu, r,
                       # priors
                       M.prior, V.prior, S.prior, df.prior,
                       # mcmc settings
                       niters, warm.up, tuning.W.tilde, tuning.phi, tuning.nu, tuning.r){
  
  
  library(fields)
  
  N.obs <- nrow(Y.obs)
  
  distmat.obs <- rdist(obs.loc)
  diameter <-  max(distmat.obs)
  
  # Initialization of K.tilde and K.tilde.inv
  
  cormat.details <- cormat.update(distmat.obs, phi, nu, r, N.obs)
  K.tilde.obs <- cormat.details$cormat.obs
  K.tilde.obs.inv <- cormat.details$cormat.obs.inv
  
  chol.V.prior <- chol(V.prior)
  prec.V.prior <- chol2inv(chol.V.prior)
  
  
  # acceptance rate for each metropolis update parameters
  
  acc.W.tilde <- acc.phi <- acc.nu <- acc.r <- 0
  
  
  W.tilde.obs.chain <- replicate(niters, matrix(NA, N.obs, q), simplify = F)
  beta.chain <- replicate(niters, matrix(NA, p, q), simplify = F)
  Sigma.chain <- replicate(niters, matrix(NA, q, q), simplify = F)
  phi.chain <- rep(NA, niters)
  nu.chain <-  rep(NA, niters)
  r.chain <-  rep(NA, niters)
  
  var.W.tilde <- matrix(1, nrow = N.obs, ncol = q)
  
  start_time <- Sys.time()  # Start time using Sys.time()
  
  for(iter in 1:niters){
    
     
    if(iter <= warm.up){
      

    W.tilde.update.details <- update.W.tilde(W.tilde.obs, beta, Sigma, K.tilde.obs,
                                             Y.obs, X.obs, N.obs, distmat.obs, 
                                             acc.W.tilde, var.W.tilde, tuning.W.tilde)
    
    W.tilde.obs <- W.tilde.update.details$W.tilde.obs
    acc.W.tilde <- W.tilde.update.details$acc.W.tilde
    
    phi.update.details <- update.phi(phi, beta, Sigma, nu, r, W.tilde.obs, distmat.obs,
                                     K.tilde.obs, K.tilde.obs.inv, acc.phi, tuning.phi)
    
    phi <- phi.update.details$phi
    K.tilde.obs <- phi.update.details$cormat.obs
    K.tilde.obs.inv <- phi.update.details$cormat.obs.inv
    acc.phi <- phi.update.details$acc.phi
  
    
    V.tilde <- chol2inv(chol(t(X.obs) %*% K.tilde.obs.inv %*% X.obs + prec.V.prior))
    M.tilde <- V.tilde %*% (t(X.obs) %*% K.tilde.obs.inv %*% W.tilde.obs + prec.V.prior %*% M.prior)
    
    chol.V.tilde <- chol(V.tilde)
    prec.V.tilde <- chol2inv(chol.V.tilde)
    
    Sigma <- update.Sigma(W.tilde.obs, X.obs, K.tilde.obs.inv, 
                          M.prior, prec.V.prior, M.tilde, prec.V.tilde, 
                          S.prior, df.prior, N.obs)
    
    beta <- update.beta(W.tilde.obs, Sigma, X.obs, K.tilde.obs.inv, M.tilde, chol.V.tilde)
    
    
    W.tilde.obs.chain[[iter]] <- W.tilde.obs
    phi.chain[iter] <- phi
    nu.chain[iter] <- nu
    r.chain[iter] <- r
    Sigma.chain[[iter]] <- Sigma
    beta.chain[[iter]] <- beta
    
    
    if(iter %% (warm.up/10) == 0) {
      elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      cat(sprintf("Progress: %.0f%%, Elapsed time: %.2f minutes\n",
                  100 * (iter / warm.up), elapsed_time))
    }
    
    
  }
  
  # After warm-up: Compute variances and adjust the tuning parameters
  if(iter == warm.up) {
    
    
    print(paste("Acceptance rate for W.tilde up to warm.up phase: ", acc.W.tilde / warm.up))
    print(paste("Acceptance rate for phi to warm.up phase: ", acc.phi / warm.up))
    
    samples <- W.tilde.obs.chain[1:warm.up]
    var.W.tilde <- matrix(NA, nrow = N.obs, ncol = q)
    
    for(i in 1:N.obs) {
      for(j in 1:q) {
        temp <- sapply(1:warm.up, function(iter) samples[[iter]][i, j])
        var.W.tilde[i, j] <- var(temp)
      }
    }
    
    var.W.tilde = var.W.tilde
    
    # Reset acceptance rates after adaptation
    acc.W.tilde <- 0
    acc.phi <- 0
    
    tuning.W.tilde <- 1.2e-2
    
    start_time <- Sys.time()  # Reset start time after adaptation

   }
  
  
    
    # Update W.tilde using the covariances after adaptation
    if (iter > warm.up) {
      
      
      W.tilde.update.details <- update.W.tilde(W.tilde.obs, beta, Sigma, K.tilde.obs,
                                               Y.obs, X.obs, N.obs, distmat.obs, 
                                               acc.W.tilde, var.W.tilde, tuning.W.tilde)
      
      W.tilde.obs <- W.tilde.update.details$W.tilde.obs
      acc.W.tilde <- W.tilde.update.details$acc.W.tilde
      
      phi.update.details <- update.phi(phi, beta, Sigma, nu, r, W.tilde.obs, distmat.obs,
                                       K.tilde.obs, K.tilde.obs.inv, acc.phi, tuning.phi)
      
      phi <- phi.update.details$phi
      K.tilde.obs <- phi.update.details$cormat.obs
      K.tilde.obs.inv <- phi.update.details$cormat.obs.inv
      acc.phi <- phi.update.details$acc.phi
      
      
      
      V.tilde <- chol2inv(chol(t(X.obs) %*% K.tilde.obs.inv %*% X.obs + prec.V.prior))
      M.tilde <- V.tilde %*% (t(X.obs) %*% K.tilde.obs.inv %*% W.tilde.obs + prec.V.prior %*% M.prior)
      
      chol.V.tilde <- chol(V.tilde)
      prec.V.tilde <- chol2inv(chol.V.tilde)
      
      Sigma <- update.Sigma(W.tilde.obs, X.obs, K.tilde.obs.inv, 
                            M.prior, prec.V.prior, M.tilde, prec.V.tilde, 
                            S.prior, df.prior, N.obs)
      
      beta <- update.beta(W.tilde.obs, Sigma, X.obs, K.tilde.obs.inv, M.tilde, chol.V.tilde)
      
      W.tilde.obs.chain[[iter]] <- W.tilde.obs
      phi.chain[iter] <- phi
      nu.chain[iter] <- nu
      r.chain[iter] <- r
      Sigma.chain[[iter]] <- Sigma
      beta.chain[[iter]] <- beta
      
      
      
      # Print progress every 10% of the remaining iterations
      if(iter %% ((niters - warm.up) / 10) == 0 || iter == niters) {
        elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
        cat(sprintf("Main Phase Progress: %.0f%%, Elapsed time: %.2f minutes\n", 
                    100 * ((iter - warm.up + 1) / (niters - warm.up)), elapsed_time))
      }

      
    } 
    
}
# Calculate total elapsed time
# end_time <- Sys.time()
# total_elapsed_time <- print(difftime(end_time, start_time, units = "mins")[[1]])  # in minutes
  
  theta.chain <- list("W.tilde.sample" = W.tilde.obs.chain,
                      "beta.sample" = beta.chain,
                      "Sigma.sample" = Sigma.chain,
                      "phi.sample" = phi.chain,
                      "nu.sample" = nu.chain,
                      "r.sample" = r.chain)
  
  # Calculate acceptance rate
  print(paste("Acceptance.W = ", acc.W.tilde/(niters-warm.up)))
  print(paste("Acceptance.phi = ", acc.phi/(niters-warm.up)))
 
  
  return(theta.chain)
  
}
