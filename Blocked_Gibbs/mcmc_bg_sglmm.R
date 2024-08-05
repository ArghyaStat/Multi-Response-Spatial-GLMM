# mcmc

mcmc.sglmm <- function(Y.obs, X.obs, obs.loc,
                       #parameters
                       W.tilde.obs, beta, Sigma, phi, nu, r,
                       # priors
                       M.prior, V.prior, S.prior, df.prior,
                       # mcmc settings
                       niters, tuning.W.tilde, tuning.phi, tuning.nu, tuning.r){
  
  
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
  
  
  W.tilde.obs.chain <- replicate(niters, matrix(NA, N, q), simplify = F)
  beta.chain <- replicate(niters, matrix(NA, p, q), simplify = F)
  Sigma.chain <- replicate(niters, matrix(NA, q, q), simplify = F)
  phi.chain <- rep(NA, niters)
  nu.chain <-  rep(NA, niters)
  r.chain <-  rep(NA, niters)
  
  W.tilde.obs.chain[[1]] <- W.tilde.obs
  phi.chain[1] <- phi
  nu.chain[1] <- nu
  r.chain[1] <- r
  Sigma.chain[[1]] <- Sigma
  beta.chain[[1]] <- beta
  
  
  
  start_time <- proc.time()[3]
  
  for(iter in 2:niters){
    
    elapsed_time <- (proc.time()[3] - start_time)/60
    
    if(iter %% ((niters)/10) == 0) {cat(sprintf("Progress: %.0f%%, Elapsed time: %.2f minutes\n", 
                                               100 * (iter / niters), elapsed_time["elapsed"]))}
    
     W.tilde.update.details <- update.W.tilde(W.tilde.obs, beta, Sigma, K.tilde.obs,
                                              Y.obs, X.obs, N.obs, distmat.obs, 
                                              acc.W.tilde, tuning.W.tilde)

     W.tilde.obs <- W.tilde.update.details$W.tilde.obs
     acc.W.tilde <- W.tilde.update.details$acc.W.tilde


     phi.update.details <- update.phi(phi, beta, Sigma, nu, r, W.tilde.obs, distmat.obs,
                                      K.tilde.obs, K.tilde.obs.inv, acc.phi, tuning.phi)

     phi <- phi.update.details$phi
     K.tilde.obs <- phi.update.details$cormat.obs
     K.tilde.obs.inv <- phi.update.details$cormat.obs.inv
     acc.phi <- phi.update.details$acc.phi



     # nu.update.details  <- update.nu(phi, beta, Sigma, nu, r, W.tilde.obs, distmat.obs,
     #                                 K.tilde.obs, K.tilde.obs.inv, acc.nu, tuning.nu)
     # 
     # nu <- nu.update.details$nu
     # K.tilde.obs <- nu.update.details$cormat.obs
     # K.tilde.obs.inv <- nu.update.details$cormat.obs.inv
     # acc.nu <- nu.update.details$acc.nu



     # r.update.details <- update.r(phi, beta, Sigma, nu, r, W.tilde.obs, distmat.obs,
     #                                K.tilde.obs, K.tilde.obs.inv, acc.r, tuning.r)
     # 
     # r <- r.update.details$r
     # K.tilde.obs <- r.update.details$cormat.obs
     # K.tilde.obs.inv <- r.update.details$cormat.obs.inv
     # acc.r <- r.update.details$acc.r
     
     Sigma <- update.Sigma(W.tilde.obs, X.obs, K.tilde.obs.inv, 
                           M.prior, prec.V.prior, M.tilde, V.tilde, 
                           S.prior, df.prior, N.obs)
     
     beta.update.details <- update.beta(W.tilde.obs, Sigma, X.obs, 
                                        K.tilde.obs.inv, M.prior, prec.V.prior)
     
     beta <- beta.update.details$beta
     M.tilde <- beta.update.details$M.tilde
     prec.V.tilde <- beta.update.details$prec.V.tilde


     W.tilde.obs.chain[[iter]] <- W.tilde.obs
     phi.chain[iter] <- phi
     nu.chain[iter] <- nu
     r.chain[iter] <- r
     Sigma.chain[[iter]] <- Sigma
     beta.chain[[iter]] <- beta
     
     
  }
  
  
  
  theta.chain <- list("W.tilde.sample" = W.tilde.obs.chain,
                      "beta.sample" = beta.chain,
                      "Sigma.sample" = Sigma.chain,
                      "phi.sample" = phi.chain,
                      "nu.sample" = nu.chain,
                      "r.sample" = r.chain)
  
  # Calculate acceptance rate.
  
  print(paste("Acceptance.W = ", acc.W.tilde/niters))
  print(paste("Acceptance.phi = ", acc.phi/niters))
  print(paste("Acceptance.nu = ", acc.nu/niters))
  print(paste("Acceptance.r = ", acc.r/niters))
  
  
  return(theta.chain)
  
  
}
