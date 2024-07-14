# mcmc

joint.mcmc.sglmm <- function(Y, X, locations,
                    # parameters   
                    W.tilde, beta, Sigma, phi, nu, r, N,
                    # hyperparams
                    M.prior, V.prior, S.prior, df.prior,
                    # mcmc settings
                    niters, tuning.W.tilde, tuning.phi, tuning.nu, tuning.r){
  
  
  library(fields)
  N <- nrow(Y)
  
  chol.V.prior <- chol(V.prior)
  prec.V.prior <- chol2inv(chol.V.prior)
  distmat <- rdist(locations)
  
  # acceptance rate for each metropolis update parameters
  
  acc.W.tilde <- acc.phi <- acc.nu <- acc.r <- 0
  
  W.tilde.chain <- replicate(niters, matrix(0, N, q), simplify = F)
  beta.chain <- replicate(niters, matrix(0, p, q), simplify = F)
  Sigma.chain <- replicate(niters, matrix(0, q, q), simplify = F)
  phi.chain <- rep(0, niters)
  nu.chain <-  rep(0, niters)
  r.chain <-  rep(0, niters)
  
  
  
  # Run Metropolis-Hastings
  for (iter in 1:niters) {
    
    if(iter %% ((niters)/10) == 0) print(paste0(100*(iter/(niters)), "%")) 
    
    ### W.tilde update
    
    W.error <- matrix(rnorm(N*q, 0, sqrt(tuning.W.tilde)), nrow = N, ncol = q)
    
    can.W.tilde <- W.tilde + W.error
    
    log.r.W.tilde <- log.posterior(can.W.tilde, beta, Sigma, phi, nu, r, Y, X, locations) - 
                     log.posterior(W.tilde, beta, Sigma, phi, nu, r, Y, X, locations)
    
    if(log(runif(1)) < log.r.W.tilde){
      
      W.tilde <- can.W.tilde
      acc.W.tilde <- acc.W.tilde + 1
      
    }
    
    ### beta_update 
    
    K.tilde <- r * Matern(distmat, 
                          range = phi, 
                          smoothness = nu) + (1- r)*diag(N)
    
    K.tilde.inv <- chol2inv(chol(K.tilde))
    
    chol.Sigma <- chol(Sigma)
    
    V.tilde <- chol2inv(chol(t(X) %*% K.tilde.inv %*% X + prec.V.prior))
    M.tilde <- V.tilde %*% (t(X) %*% K.tilde.inv %*% W.tilde + prec.V.prior %*% M.prior)
    
    chol.V.tilde <- chol(V.tilde)
    
    beta.samp <- matrix(rnorm(p*q, 0, 1), nrow = p , ncol = q)
    
    beta <- M.tilde + t(chol.V.tilde) %*% beta.samp %*% chol.Sigma
    
    beta
    
    
    ### Sigma Update
    
    resid.W.tilde <- W.tilde - X %*% beta
    rss.W.tilde <- t(resid.W.tilde) %*% K.tilde.inv %*% resid.W.tilde
    
    rss.beta <- t(beta - M.prior) %*% prec.V.prior %*% (beta - M.prior)
    
    S.tilde <- S.prior + rss.W.tilde + rss.beta
    df.tilde <- df.prior + N + p
    
    Sigma <- riwish(v = df.tilde, S = S.tilde)
    
    Sigma
    
    ### phi update
    
    can.phi <- rnorm(1, phi, sqrt(tuning.phi))
    
    # Compute log posterior for the proposed value
    log.r.phi <- log.posterior(W.tilde, beta, Sigma, can.phi, nu, r, Y, X, locations) - 
                 log.posterior(W.tilde, beta, Sigma, phi, nu, r, Y, X, locations)
    
    # Accept or reject
    if (log(runif(1)) < log.r.phi) {
      
      phi <- can.phi
      acc.phi <- acc.phi + 1
      
    }
    
    ### nu update
    
    can.nu <- rnorm(1, nu, sqrt(tuning.nu))
    
    
    log.r.nu <- log.posterior(W.tilde, beta, Sigma, phi, can.nu, r, Y, X, locations) - 
                log.posterior(W.tilde, beta, Sigma, phi, nu, r, Y, X, locations)
    
    # Accept or reject
    if (log(runif(1)) < log.r.nu) {
      
      nu <- can.nu
      acc.nu <- acc.nu + 1
      
    }
    
    ### r update
    
    can.r <- rnorm(1, r, sqrt(tuning.r))
    
    
    log.r.r <- log.posterior(W.tilde, beta, Sigma, phi, nu, can.r, Y, X, locations) - 
               log.posterior(W.tilde, beta, Sigma, phi, nu, r, Y, X, locations)
    
    # Accept or reject
    if (log(runif(1)) < log.r.r) {
      
      r <- can.r
      acc.r <- acc.r + 1
      
    }
    
    
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
