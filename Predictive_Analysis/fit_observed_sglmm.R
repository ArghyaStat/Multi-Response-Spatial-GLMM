rm(list = ls())
#
library(this.path)
library(fields)
library(MCMCpack)
library(rlist)

mydir <- this.path::here()
setwd(mydir)

load("sglmm.simulation.Rdata")

#--------------------

source("update_parameters_sglmm.R")
source("mcmc_sglmm.R")

# Sample initial parameters at true value
W.tilde.obs <- true.W.tilde.obs
  # matrix(rep(0, N*q), nrow = N, ncol = q)
beta <- true.beta
# matrix(rep(0, p*q), nrow = p, ncol = q)
Sigma <- true.Sigma
#diag(c(3,4))
phi <- true.phi
nu <- true.nu
r <- true.r

# Number of iterations
niters <- 1e4

# Tuning parameters list
tuning.W.tilde <- 6e-3
tuning.phi <- 1e-2
tuning.nu <- 8e-2
tuning.r <- 1e-2

# Separate Tuning parameters list
# tuning.W.tilde <- 6e-3
# tuning.phi <- 1e-2
# tuning.nu <- 8e-2 
# tuning.r <- 1e-2

distmat.obs <- rdist(obs.loc)


M.prior <- matrix(0, p, q)
V.prior <- 1e2 * diag(p) 
S.prior <- diag(q)
df.prior <- q + 1

fit.observed.sglmm <- mcmc.sglmm(Y.obs, X.obs, obs.loc,
                                 W.tilde.obs, beta, Sigma, phi, nu, r,
                                 # priors
                                 M.prior, V.prior, S.prior, df.prior,
                                 # mcmc settings
                                 niters, tuning.W.tilde, tuning.phi, tuning.nu, tuning.r)

# Saving MCMC chain
list.save(fit.observed.sglmm, file = "separable.mgp.MCMC.chain.Rdata")

theta.chain <- list.load(file = "separable.mgp.MCMC.chain.Rdata")

# Traceplots

par(mfrow = c(1,1))

pdf("traceplot_phi.pdf", width=5 ,height=5)
trace.phi <- plot.ts(theta.chain$phi.sample, ylab = "phi", main = "Traceplot of phi")
abline(h = true.phi, col = 'blue', lwd = 2)
dev.off()

pdf("traceplot_nu.pdf", width=5 ,height=5)
trace.nu <- plot.ts(theta.chain$nu.sample, ylab = "nu", main = "Traceplot of nu")
abline(h = true.nu, col = 'blue', lwd = 2)
dev.off()

pdf("traceplot_r.pdf", width=5 ,height=5)
trace.r <- plot.ts(theta.chain$r.sample, ylab = "r", main = "Traceplot of r")
abline(h = true.r, col = 'blue', lwd = 2)
dev.off()

# acfplots

pdf("acf_phi.pdf", width=5 ,height=5)
acf.phi <- acf(theta.chain$phi.sample, main = "ACF plot of phi", lag.max = 200)
dev.off()

pdf("acf_nu.pdf", width=5 ,height=5)
acf.nu <- acf(theta.chain$nu.sample, main = "ACF plot of nu", lag.max = 200)
dev.off()

pdf("acf_r.pdf", width=5 ,height=5)
acf.r <- acf(theta.chain$r.sample, main = "ACF plot of r", lag.max = 200)
dev.off()




labels.beta <- list()
for (i in 1:p) {
  for (j in 1:q) {
    labels.beta <- c(labels.beta, 
                     list(paste0('beta (', i, ',', j, ')')))
  }
}


pdf("traceplot_beta.pdf", width=7.5 ,height=6)
par(mfrow=c(p,q))

for (i in 1:p) {
  for (j in 1:q) {
    plot(1:niters, sapply(theta.chain$beta.sample, function(x) x[i, j]), 
         type='l', col=1, xlab='Iteration', ylab='beta',
         main=labels.beta[(i-1)*q + j])
    abline(h = true.beta[i,j], col = 'blue', lwd = 2)
    
  }
}
dev.off()

pdf("acf_beta.pdf", width= 7.5 ,height=6)
par(mfrow=c(p,q))

for (i in 1:p) {
  for (j in 1:q) {
    
    acf(sapply(theta.chain$beta.sample, function(x) x[i, j]), 
        main = labels.beta[(i-1)*q + j], lag.max = 100)
    
  }
}
dev.off()

# Output Analysis of Sigma

labels.Sigma <- list()
for (i in 1:p) {
  for (j in 1:q) {
    labels.Sigma <- c(labels.Sigma, 
                      list(paste0('Sigma (', i, ',', j, ')')))
  }
}


pdf("traceplot_Sigma.pdf", width = 7.5 ,height = 6)
par(mfrow=c(q,q))

for (i in 1:q) {
  for (j in 1:q) {
    plot(1:niters, sapply(theta.chain$Sigma.sample, function(x) x[i, j]), 
         type='l', col=1, xlab='Iteration', ylab='Sigma',
         main=labels.Sigma[(i-1)*q + j])
    abline(h = true.Sigma[i,j], col = 'blue', lwd = 2)
    
  }
}
dev.off()


pdf("acf_Sigma.pdf", width = 7.5 ,height = 6)
par(mfrow=c(q,q))

for (i in 1:q) {
  for (j in 1:q) {
    
    acf(sapply(theta.chain$Sigma.sample, function(x) x[i, j]), 
        main=labels.Sigma[(i-1)*q + j],
        lag.max = 200)
    
  }
}
dev.off()

