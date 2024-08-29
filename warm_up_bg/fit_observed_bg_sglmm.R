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
  niters <- 1.1e5
  warm.up <- 1e4
  
  
  # Tuning parameters list for log(1+BA) and CNT
  tuning.W.tilde <- 3.6e-3
  tuning.phi <- 2.8e-3
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
list.save(theta.chain, file = "bg.sglmm.MCMC.chain.Rdata")

theta.chain <- list.load(file = "bg.sglmm.MCMC.chain.Rdata")



# variance computation for w.tilde till warm up

samples <- theta.chain$W.tilde.sample[1:warm.up]
var.W.tilde <- matrix(NA, nrow = N.obs, ncol = q)

for(i in 1:N.obs) {
  for(j in 1:q) {
    temp <- sapply(1:warm.up, function(iter) samples[[iter]][i, j])
    var.W.tilde[i, j] <- var(temp)
  }
}
summary(var.W.tilde)



#set.seed(100)  # For reproducibility
selected_rows <- sample(1:N.obs, 5, replace = FALSE)


# Output Analysis of W.tilde

pdf("traceplot_W.tilde.pdf", width = 6*q, height = 4*q)
par(mfrow=c(p, q), mar = c(5,5,3,1))

for (i in selected_rows) {
  for (j in 1:q) {
    plot(1:(niters - warm.up), xlim = c(1, (niters - warm.up)), 
         sapply(theta.chain$W.tilde.sample[(warm.up + 1) : niters], function(x) x[i, j]), 
         type='l', col=1, xlab='Iterations', 
         ylab=bquote(W.tilde[.(i) * "," * .(j)]),
         main=bquote("Traceplot of " * W.tilde[.(i) * "," * .(j)]))
    abline(h = W.tilde.obs[i,j], col = 'blue', lwd = 2)
  }
}

dev.off()

pdf("acf_W.tilde.pdf", width = 6*q, height = 4*q)
par(mfrow=c(p, q), mar = c(5,5,3,1))

for (i in 1:selected_rows) {
  for (j in 1:q) {
    
    acf(sapply(theta.chain$W.tilde.sample[(warm.up + 1) : niters], function(x) x[i, j]), 
        ylab=bquote(W.tilde[.(i) * "," * .(j)]),
        main=bquote("ACF plot of " * W.tilde[.(i) * "," * .(j)]),
        lag.max = 500)
    
  }
}
dev.off()




# Traceplots

par(mfrow = c(1,1))

pdf("traceplot_phi.pdf", width=4 ,height=4)
trace.phi <- plot.ts(theta.chain$phi.sample[(warm.up + 1) : niters], ylab = expression(phi), xlab = "Iterations",
                     main = expression("Traceplot of " * phi))
abline(h = true.phi, col = 'blue', lwd = 2)
dev.off()



# acfplots

pdf("acf_phi.pdf", width=4 ,height=4)
acf.phi <- acf(theta.chain$phi.sample[(warm.up + 1) : niters], main = expression("ACF plot of " * phi),  
               ylab = expression(phi), lag.max = 100)
dev.off()








# Output Analysis of beta

pdf("traceplot_beta.pdf", width = 6*q, height = 4*q)
par(mfrow=c(p, q), mar = c(5,5,3,1))

for (i in 1:p) {
  for (j in 1:q) {
    plot(1:(niters - warm.up), xlim = c(1, (niters - warm.up)), 
         sapply(theta.chain$beta.sample[(warm.up + 1) : niters], function(x) x[i, j]), 
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
    
    acf(sapply(theta.chain$beta.sample[(warm.up + 1) : niters], function(x) x[i, j]), 
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
    plot(1:(niters - warm.up), xlim = c(1, (niters-warm.up)), 
         sapply(theta.chain$Sigma.sample[(warm.up + 1) : niters], function(x) x[i, j]), 
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
    
    acf(sapply(theta.chain$Sigma.sample[(warm.up + 1) : niters], function(x) x[i, j]), 
        ylab=bquote(Sigma[.(i) * "," * .(j)]),
        main=bquote("ACF plot of " * Sigma[.(i) * "," * .(j)]),
        lag.max = 100)
    
  }
}
dev.off()

# pdf("traceplot_nu.pdf", width=4 ,height=4)
# trace.nu <- plot.ts(theta.chain$nu.sample, ylab = expression(nu), xlab = "Iterations",
#                     main = expression("Traceplot of " * nu))
# abline(h = true.nu, col = 'blue', lwd = 2)
# dev.off()
# 
# pdf("traceplot_r.pdf", width=4 ,height=4)
# trace.r <- plot.ts(theta.chain$r.sample, ylab = expression(r), xlab = "Iterations",
#                    main = expression("Traceplot of " * r))
# abline(h = true.r, col = 'blue', lwd = 2)
# dev.off()
# 
# 
# pdf("acf_nu.pdf", width=4 ,height=4)
# acf.nu <- acf(theta.chain$nu.sample, main = expression("ACF plot of " * nu), 
#               ylab = expression(nu),
#               lag.max = 100)
# dev.off()
# 
# pdf("acf_r.pdf", width=4 ,height=4)
# acf.r <- acf(theta.chain$r.sample, main = expression("ACF plot of " * r), 
#              ylab = expression(r),
#              lag.max = 100)
# dev.off()