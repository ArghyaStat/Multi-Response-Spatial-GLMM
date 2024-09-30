  rm(list = ls())
  #
  library(this.path)
  library(fields)
  library(fBasics)
  library(MCMCpack)
  library(truncnorm)
  #library(mvtnorm)  # one can cross-verify dmvn with dmvnorm
  library(rlist)
  
  mydir <- this.path::here()
  setwd(mydir)
  
  load("sglmm.bg.simulation.Rdata")
  
  #--------------------
  
  source("update_parameters_bg_appx_sglmm.R")
  source("mcmc_bg_appx_sglmm.R")
  
  
  # Number of iterations
  niters <- 1e4
  warm.up <- 0  ### no warming up here
  
  
  # Tuning parameters list
  
  tuning.W.tilde <- 5e-2   # chosen from 1e-20 to 5.
  tuning.W.tilde <- 1e-2
  tuning.phi <- 1e-3
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
list.save(theta.chain, file = "sglmm.bg.appx.chain.Rdata")

theta.chain <- list.load(file = "sglmm.bg.appx.chain.Rdata")

## randomly select 5 sites out of N.obs for traceplots

#set.seed(100)  # For reproducibility
selected_rows <- sample(1:N.obs, 5, replace = FALSE)


# Output Analysis of W.tilde

pdf("traceplot_W.tilde.pdf", width = 6*q, height = 4*q)
par(mfrow=c(p, q), mar = c(5,5,3,1))


# tarceplot of W
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

# Function to compute the ACF for a matrix
compute_acf_matrix <- function(chain, i, j, niters) {
  component_values <- numeric(niters)
  
  for (k in 1:niters) {
    component_values[k] <- chain[[k]][i, j]
  }
  
  acf(component_values, plot = FALSE, lag.max = 500)
}

# Function to compute the ACF for a numeric value
compute_acf_numeric <- function(chain, niters) {
  values <- numeric(niters)
  
  for (k in 1:niters) {
    values[k] <- chain[[k]]
  }
  
  acf(values, plot = FALSE, lag.max = 500)
}

acf.W_tilde <- compute_acf_matrix(theta.chain$W.tilde.sample[(warm.up + 1) : niters], i, j, (niters - warm.up))


#### acf for W

pdf("acf_W.tilde.pdf", width = 6*q, height = 4*q)
par(mfrow=c(p, q), mar = c(5,5,3,1))

for (i in selected_rows) {
  for (j in 1:q) {
   
    plot(acf.W_tilde$lag, acf.W_tilde$acf, type = "l", lwd = 2, col = "black",
         ylab = bquote(W[tilde][.(i) * "," * .(j)]),
         main = bquote("ACF plot of " * W.tilde[.(i) * "," * .(j)]),
         xlab = "Lag", ylim = range(c(0,1)))
    
  }
}
dev.off()

# Traceplot for phi

par(mfrow = c(1,1))

pdf("traceplot_phi.pdf", width=4 ,height=4)
trace.phi <- plot.ts(theta.chain$phi.sample[(warm.up + 1) : niters], ylab = expression(phi), xlab = "Iterations",
                     main = expression("Traceplot of " * phi))
abline(h = true.phi, col = 'blue', lwd = 2)
dev.off()


### acf computation for phi

acf.phi <- compute_acf_numeric(theta.chain$phi, (niters - warm.up))


# acfplot for phi

pdf("acf_phi.pdf", width=4 ,height=4)

plot(acf.phi$lag, acf.phi$acf, type = "l", lwd = 2, col = "black",
     ylab = expression(phi),
     main = expression("ACF plot of " * phi),
     xlab = "Lag", ylim = c(0,1))


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

# acf plot for beta

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

# acf plot for Sigma

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


# post.warmup.bg.chain <- list("W.tilde.sample" = theta.chain$W.tilde.sample[(warm.up + 1): niters],
#                              "beta.sample" = theta.chain$beta.sample[(warm.up + 1): niters],
#                              "Sigma.sample" = theta.chain$Sigma.sample[(warm.up + 1): niters],
#                              "phi.sample" = theta.chain$phi.sample[(warm.up + 1): niters]
# )
# 
# 
# 
# 
# list.save(post.warmup.bg.chain, file = "post.warm.up.bg.chain.Rdata")
# 
# abc <- list.load(file =  "post.warm.up.bg.chain.Rdata")
# 
# # variance computation for w.tilde till warm up
# 
# samples <- theta.chain$W.tilde.sample[1:warm.up]
# var.W.tilde <- matrix(NA, nrow = N.obs, ncol = q)
# 
# for(i in 1:N.obs) {
#   for(j in 1:q) {
#     temp <- sapply(1:warm.up, function(iter) samples[[iter]][i, j])
#     var.W.tilde[i, j] <- var(temp)
#   }
# }
# summary(var.W.tilde)
