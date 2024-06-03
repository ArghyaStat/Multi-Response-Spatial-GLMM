# mcmc

mcmc.sglmm <- function(Y, X, locations,
                    W.tilde, beta, Sigma, phi, nu, r, N,
                    # priors
                    M.prior, V.prior, S.prior, df.prior, distmat,
                    # mcmc settings
                    niters, tuning.W.tilde, tuning.phi, tuning.nu, tuning.r){
  
  
  library(fields)
  
  
  N <- nrow(Y)
  
  d <- rdist(locations)
  
  # Initialization of K.tilde and K.tilde.inv
  
  cormat.details <- cormat.update(d, phi, nu, r, N)
  K.tilde <- cormat.details$cormat
  K.tilde.inv <- cormat.details$cormat.inv
  
  # acceptance rate for each metropolis update parameters
  
  acc.W.tilde <- acc.phi <- acc.nu <- acc.r <- 0
  
  W.tilde.chain <- replicate(niters, matrix(0, N, q), simplify = F)
  beta.chain <- replicate(niters, matrix(0, p, q), simplify = F)
  Sigma.chain <- replicate(niters, matrix(0, q, q), simplify = F)
  phi.chain <- rep(0, niters)
  nu.chain <-  rep(0, niters)
  r.chain <-  rep(0, niters)
  
  beta <- beta
  Sigma <- Sigma 
  phi <- phi
  nu <- nu
  r <- r
  
  for(iter in 1:niters){
    
    if(iter %% ((niters)/10) == 0) print(paste0(100*(iter/(niters)), "%"))
     
    
     W.tilde.update.details <- update.W.tilde(beta, Sigma, r, phi, nu, Y, X, N, K.tilde,
                                             distmat, acc.W.tilde, tuning.W.tilde)
    
     W.tilde <- W.tilde.update.details$W.tilde
     acc.W.tilde <- W.tilde.update.details$acc.W.tilde
     
     #print(W.tilde)
     
     Sigma <- update.Sigma(W.tilde, X, K.tilde.inv,
                           M.prior, V.prior, S.prior, df.prior, N)
     
     #print(Sigma)

     beta <- update.beta(W.tilde, Sigma, X, K.tilde.inv, M.prior, V.prior)
    
     #print(beta)
    
     phi.update.details <- update.phi(phi, beta, Sigma, nu, r, W.tilde, distmat,
                                      K.tilde, K.tilde.inv, acc.phi, tuning.phi)

     phi <- phi.update.details$phi
     K.tilde <- phi.update.details$cormat
     K.tilde.inv <- phi.update.details$cormat.inv
     acc.phi <- phi.update.details$acc.phi
     
     #print(phi)

     nu.update.details  <- update.nu(phi, beta, Sigma, nu, r, W.tilde, distmat,
                                     K.tilde, K.tilde.inv, acc.nu, tuning.nu)

     nu <- nu.update.details$nu
     K.tilde <- nu.update.details$cormat
     K.tilde.inv <- nu.update.details$cormat.inv
     acc.nu <- nu.update.details$acc.nu
     
     #print(nu)

     r.update.details <- update.r(phi, beta, Sigma, nu, r, W.tilde, distmat,
                                    K.tilde, K.tilde.inv, acc.r, tuning.r)

     r <- r.update.details$r
     K.tilde <- r.update.details$cormat
     K.tilde.inv <- r.update.details$cormat.inv
     acc.r <- r.update.details$acc.r
     
     #print(r)

    

     W.tilde.chain[[iter]] <- W.tilde
     beta.chain[[iter]] <- beta
     Sigma.chain[[iter]] <- Sigma
     phi.chain[iter] <- phi
     nu.chain[iter] <- nu
     r.chain[iter] <- r
     
  }
  
  theta.chain <- list("W.tilde.sample" = W.tilde.chain,
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
