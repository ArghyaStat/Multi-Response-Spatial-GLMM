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
W.tilde <- W.mat
  # matrix(rep(0, N*q), nrow = N, ncol = q)
beta <- true.beta
# matrix(rep(0, p*q), nrow = p, ncol = q)
Sigma <- true.Sigma
phi <- true.phi
nu <- true.nu
r <- true.r

# Number of iterations
niters <- 1e4

# Tuning parameters list

tuning.W.tilde <- 5e-3
tuning.phi <- 1e-2
tuning.nu <- 8e-2 
tuning.r <- 5e-2

distmat <- rdist(locations)

p <- 3
q <- 2
N <- 5e2

M.prior <- matrix(0, p, q)
V.prior <- 1e2 * diag(p) 
S.prior <- diag(q)
df.prior <- q + 1

fit.observed.sglmm <- mcmc.sglmm(Y, X, locations,
                           W.tilde, beta, Sigma, phi, nu, r, N,
                           # priors
                           M.prior, V.prior, S.prior, df.prior,
                           # mcmc settings
                           niters, tuning.W.tilde, tuning.phi, tuning.nu, tuning.r)

# Saving MCMC chain
list.save(fit.observed.sglmm, file = "separable.mgp.MCMC.chain.Rdata")

theta.chain <- list.load(file = "separable.mgp.MCMC.chain.Rdata")

# Traceplots

par(mfrow = c(1,1))
trace.phi <- plot.ts(theta.chain$phi.sample, ylab = "phi", main = "Traceplot of phi")
abline(h = true.phi, col = 'blue', lwd = 2)

trace.nu <- plot.ts(theta.chain$nu.sample, ylab = "nu", main = "Traceplot of nu")
abline(h = true.nu, col = 'blue', lwd = 2)

trace.r <- plot.ts(theta.chain$r.sample, ylab = "r", main = "Traceplot of r")
abline(h = true.r, col = 'blue', lwd = 2)

# acfplots

acf.phi <- acf(theta.chain$phi.sample, main = "ACF plot of phi", lag.max = 200)
acf.nu <- acf(theta.chain$nu.sample, main = "ACF plot of nu", lag.max = 200)
acf.r <- acf(theta.chain$r.sample, main = "ACF plot of r", lag.max = 200)





labels.beta <- list()
for (i in 1:p) {
  for (j in 1:q) {
    labels.beta <- c(labels.beta, 
                     list(paste0('beta (', i, ',', j, ')')))
  }
}



par(mfrow=c(p,q))

for (i in 1:p) {
  for (j in 1:q) {
    plot(1:niters, sapply(theta.chain$beta.sample, function(x) x[i, j]), 
         type='l', col=1, xlab='Iteration', ylab='beta',
         main=labels.beta[(i-1)*q + j])
    abline(h = true.beta[i,j], col = 'blue', lwd = 2)
    
  }
}

par(mfrow=c(p,q))

for (i in 1:p) {
  for (j in 1:q) {
    
    acf(sapply(theta.chain$beta.sample, function(x) x[i, j]), 
        main = labels.beta[(i-1)*q + j], lag.max = 100)
    
  }
}


# Output Analysis of Sigma

labels.Sigma <- list()
for (i in 1:p) {
  for (j in 1:q) {
    labels.Sigma <- c(labels.Sigma, 
                      list(paste0('Sigma (', i, ',', j, ')')))
  }
}

par(mfrow=c(q,q))

for (i in 1:q) {
  for (j in 1:q) {
    plot(1:niters, sapply(theta.chain$Sigma.sample, function(x) x[i, j]), 
         type='l', col=1, xlab='Iteration', ylab='Sigma',
         main=labels.Sigma[(i-1)*q + j])
    abline(h = true.Sigma[i,j], col = 'blue', lwd = 2)
    
  }
}

par(mfrow=c(q,q))

for (i in 1:q) {
  for (j in 1:q) {
    
    acf(sapply(theta.chain$Sigma.sample, function(x) x[i, j]), 
        main=labels.Sigma[(i-1)*q + j],
        lag.max = 200)
    
  }
}
