rm(list = ls())
#
library(this.path)
library(fields)
library(MCMCpack)
library(rlist)

mydir <- this.path::here()
setwd(mydir)

load("wildfire.Rdata")

#--------------------

source("update_parameters_data.R")
source("mcmc_data.R")


# Number of iterations
niters <- 1e5

# Tuning parameters list univariate BA
# tuning.W.tilde <- 2.5e-3
# tuning.phi <- 4e-1
# tuning.nu <- 2e-2
# tuning.r <- 1e-2

# Tuning parameters list univariate CNT
tuning.W.tilde <- 3e-3
tuning.phi <- 4.8e-1
tuning.nu <- 2e-2
tuning.r <- 1e-2

# Tuning parameters list for log(1+BA) and CNT
# tuning.W.tilde <- 3e-3
# tuning.phi <- 0.12
# tuning.nu <- 2e-2
# tuning.r <- 1e-2

# Tuning parameters list for log(1+BA) and log(1+CNT)
# tuning.W.tilde <- 1.5e-4
# tuning.phi <- 0.8
# tuning.nu <- 2e-2
# tuning.r <- 1e-2

distmat.obs <- rdist(obs.loc)

M.prior <- matrix(0, p, q)
V.prior <-  1e2*diag(p) 
S.prior <-  diag(q)
df.prior <- q + 1

# Sample initial parameters at true value
phi <- 0.1 * diameter
nu <- 0.5
r <- 1
beta <- (chol2inv(chol((t(X.obs) %*% X.obs)))) %*% (t(X.obs) %*% Y.obs)
W.tilde.obs <- X.obs %*% beta
Sigma <- diag(q)
 

fit.observed.data <- mcmc.sglmm(Y.obs, X.obs, obs.loc,
                                 W.tilde.obs, beta, Sigma, phi, nu, r,
                                 # priors
                                 M.prior, V.prior, S.prior, df.prior,
                                 # mcmc settings
                                 niters, tuning.W.tilde, tuning.phi, tuning.nu, tuning.r)

# Saving MCMC chain
list.save(fit.observed.data, file = "wildfire.MCMC.chain.Rdata")

theta.chain <- list.load(file = "wildfire.MCMC.chain.Rdata")




# Traceplots

par(mfrow = c(1,1))

pdf("traceplot_phi.pdf", width=5 ,height=5)
trace.phi <- plot.ts(theta.chain$phi.sample, ylab = "phi", main = "Traceplot of phi")
dev.off()

pdf("traceplot_nu.pdf", width=5 ,height=5)
trace.nu <- plot.ts(theta.chain$nu.sample, ylab = "nu", main = "Traceplot of nu")
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
for (i in 1:q) {
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

