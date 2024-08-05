  rm(list = ls())
  #
  library(this.path)
  library(fields)
  library(MCMCpack)
  library(rlist)
  
  mydir <- this.path::here()
  setwd(mydir)
  
  load("sglmm.bg.simulation.Rdata")
  
  #--------------------
  
  source("update_parameters_bg_sglmm.R")
  source("mcmc_bg_sglmm.R")
  
  
  # Number of iterations
  niters <- 1e5
  
  
  # Tuning parameters list for log(1+BA) and CNT
  tuning.W.tilde <- 3.5e-3
  tuning.phi <- 3e-3
  tuning.nu <- 2e-2
  tuning.r <- 1e-2
  
  
  distmat.obs <- rdist(obs.loc)
  
  M.prior <- matrix(0, p, q)
  V.prior <-  1e2*diag(p) 
  S.prior <-  diag(q)
  df.prior <- q + 1
  
  # Sample initial parameters at true value
  phi <- 0.1 * diameter
  nu <- 0.5
  r <- 1
  W.tilde.obs <- cbind(Y.obs[,1], log(1+ Y.obs[,2]))
  beta <- (chol2inv(chol((t(X.obs) %*% X.obs)))) %*% (t(X.obs) %*% W.tilde.obs)
  Sigma <- t(W.tilde.obs - X.obs %*% beta) %*% (W.tilde.obs - X.obs %*% beta)/N.obs
  
   
  
  fit.observed.data <- mcmc.sglmm(Y.obs, X.obs, obs.loc,
                                   W.tilde.obs, beta, Sigma, phi, nu, r,
                                   # priors
                                   M.prior, V.prior, S.prior, df.prior,
                                   # mcmc settings
                                   niters, tuning.W.tilde, tuning.phi, tuning.nu, tuning.r)

# Saving MCMC chain
list.save(fit.observed.data, file = "bg.sglmm.MCMC.chain.Rdata")

theta.chain <- list.load(file = "bg.sglmm.MCMC.chain.Rdata")




# Traceplots

par(mfrow = c(1,1))

pdf("traceplot_phi.pdf", width=4 ,height=4)
trace.phi <- plot.ts(theta.chain$phi.sample, ylab = expression(phi), xlab = "Iterations",
                     main = expression("Traceplot of " * phi))
abline(h = true.phi, col = 'blue', lwd = 2)
dev.off()

pdf("traceplot_nu.pdf", width=4 ,height=4)
trace.nu <- plot.ts(theta.chain$nu.sample, ylab = expression(nu), xlab = "Iterations",
                    main = expression("Traceplot of " * nu))
abline(h = true.nu, col = 'blue', lwd = 2)
dev.off()

pdf("traceplot_r.pdf", width=4 ,height=4)
trace.r <- plot.ts(theta.chain$r.sample, ylab = expression(r), xlab = "Iterations",
                   main = expression("Traceplot of " * r))
abline(h = true.r, col = 'blue', lwd = 2)
dev.off()

# acfplots

pdf("acf_phi.pdf", width=4 ,height=4)
acf.phi <- acf(theta.chain$phi.sample, main = expression("ACF plot of " * phi),  
               ylab = expression(phi), lag.max = 100)
dev.off()

pdf("acf_nu.pdf", width=4 ,height=4)
acf.nu <- acf(theta.chain$nu.sample, main = expression("ACF plot of " * nu), 
              ylab = expression(nu),
              lag.max = 100)
dev.off()

pdf("acf_r.pdf", width=4 ,height=4)
acf.r <- acf(theta.chain$r.sample, main = expression("ACF plot of " * r), 
             ylab = expression(r),
             lag.max = 100)
dev.off()



# Output Analysis of beta

pdf("traceplot_beta.pdf", width = 6*q, height = 4*q)
par(mfrow=c(p, q), mar = c(5,5,3,1))

for (i in 1:p) {
  for (j in 1:q) {
    plot(1:niters, xlim = c(1, niters), 
         sapply(theta.chain$beta.sample, function(x) x[i, j]), 
         type='l', col=1, xlab='Iterations', 
         ylab=bquote(beta[.(i) * "," * .(j)]),
         main=bquote("Traceplot of " * beta[.(i) * "," * .(j)]))
    abline(h = true.beta[i,j], col = 'blue', lwd = 2)
  }
}

dev.off()


pdf("acf_beta.pdf", width = 6*q, height = 4*q)
par(mfrow=c(p, q), mar = c(5,5,3,1))

for (i in 1:p) {
  for (j in 1:q) {
    
    acf(sapply(theta.chain$beta.sample, function(x) x[i, j]), 
        ylab=bquote(beta[.(i) * "," * .(j)]),
        main=bquote("ACF plot of " * beta[.(i) * "," * .(j)]),
        lag.max = 100)
    
  }
}
dev.off()

# Output Analysis of Sigma

pdf("traceplot_Sigma.pdf", width = 6*q, height = 4*q)
par(mfrow=c(q, q))

for (i in 1:q) {
  for (j in 1:q) {
    plot(1:niters, xlim = c(1, niters), 
         sapply(theta.chain$Sigma.sample, function(x) x[i, j]), 
         type='l', col=1, xlab='Iterations', 
         ylab=bquote(Sigma[.(i) * "," * .(j)]),
         main=bquote("Traceplot of " * Sigma[.(i) * "," * .(j)]))
    abline(h = true.Sigma[i,j], col = 'blue', lwd = 2)
  }
}

dev.off()



pdf("acf_Sigma.pdf", width = 6*q ,height = 4*q)
par(mfrow=c(q,q))

for (i in 1:q) {
  for (j in 1:q) {
    
    acf(sapply(theta.chain$Sigma.sample, function(x) x[i, j]), 
        ylab=bquote(Sigma[.(i) * "," * .(j)]),
        main=bquote("ACF plot of " * Sigma[.(i) * "," * .(j)]),
        lag.max = 100)
    
  }
}
dev.off()

