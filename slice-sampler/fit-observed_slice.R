rm(list = ls())

library(spam)
library(this.path)
library(fields)
library(fBasics)
library(MCMCpack)
library(truncnorm)
library(mvtnorm)  # one can cross-verify dmvn with dmvnorm
library(rlist)
library(foreach)
library(doParallel)
library(scoringRules)
library(FNN)
library(mcmcse)

libraries <- c("fields", "fBasics", "MCMCpack", "truncnorm", 
               "rlist", "foreach","doParallel", "spam")
functions <- c("dmatnorm.sgv", "dmvn.sgv", "log.likelihood", "update.W", 
               "update.beta", "update.Sigma", "update.phi", "max_min",
               "neighbor_matrix", "comb", "dist.nn", "U.sgv","L_inv", 
               "tuning.update")
mydir <- this.path::here()
setwd(mydir)

source("data_simulation.R")
source("vecchia_slice.R")
source("aux_functions_slice.R")
source("update_parameters_slice.R")
source("warmup_slice.R")
source("mcmc_main_slice.R")
source("pred_vecchia_slice.R")
source("predictive_scores_slice.R")

set.seed(1709)
p <- 3
q <- 2
true.beta <- matrix(rnorm(p * q, mean = 0, sd = 1), nrow = p, ncol = q)
true.Sigma <- matrix(c(2,0.8,0.8,1), nrow = q, ncol = q, byrow = TRUE)
true.phi <- 0.1
true.nu <- 0.5
pred.prop <- 0.2

data <- sim.data(q = 2, N = 2e2, 
                 family = c("Gaussian", "Poisson"),
                 true.beta = true.beta,
                 true.Sigma = true.Sigma, 
                 true.phi = true.phi,
                 true.nu = true.nu,
                 pred.prop = 0.2)

family <- c("Gaussian", "Poisson")

#### Prior specifications ####

M.prior <- matrix(0, p, q)
V.prior <-  1e2*diag(p) 
S.prior <-  diag(q)
df.prior <- q + 1


### Vecchia specifications and pre-computations ###

m <- 20
obs.locs <- data$obs.locs
N.obs <- data$N.obs
Y.obs <- data$Y.obs
X.obs <- data$X.obs


obs.ord <- max_min(obs.locs) 
 # Assume max_min returns max-min order indices
obs.locs.ord <- obs.locs[obs.ord, , drop = FALSE]  # Reorder locations based on max-min ordering
distobs.ord <- rdist(obs.locs.ord)
diameter <-  max(distobs.ord)
b_phi <- diameter/log(10)


NNarray.obs <- neighbor_matrix(obs.locs.ord, m)
distobs.nn <- dist.nn(obs.locs.ord, neighbor_matrix = NNarray.obs)

# Reorder Y.obs, X.obs, and W.obs accordingly
Y.obs.ord <- Y.obs[obs.ord, , drop = FALSE]
X.obs.ord <- X.obs[obs.ord, , drop = FALSE]

# Sample initial parameters at true value
phi <- data$true.phi
nu <- data$true.nu
W.obs <- cbind(Y.obs[,1], log(1 + Y.obs[,2]))
W.obs.ord <- W.obs[obs.ord, , drop = FALSE]
beta <- chol2inv(chol(crossprod(X.obs.ord))) %*% crossprod(X.obs.ord, W.obs.ord)
chol.K.obs.ord.inv <- U.sgv(distobs.nn, NNarray.obs, 
                              phi, nu, m)
rss.W.part <- chol.K.obs.ord.inv %*% (W.obs - X.obs %*% beta)
Sigma <- crossprod.spam(rss.W.part) /N.obs

Sigma.inv <- chol2inv(chol(Sigma))
chol.Sigma.inv <- chol(Sigma.inv)
chol.Sigma <- chol(Sigma)

tuning.phi <- 3e-3

# Number of iterations
niters <- 1e4



fit.main <- mcmc.main(Y.obs.ord, X.obs.ord, obs.locs.ord, m, N.obs, p, q,
                      distobs.ord, distobs.nn, NNarray.obs, 
                      family, nu, W.obs.ord, beta, Sigma, phi,
                      M.prior, V.prior, S.prior, df.prior, b_phi,
                      niters, tuning.phi)

burnin <- 0
pred.iters <- niters

W.ord.post <- fit.main$W.ord.samples
beta.post <- fit.main$beta.samples
Sigma.post <- fit.main$Sigma.samples
chol.Sigma.post <- fit.main$chol.Sigma.samples
phi.post <- fit.main$phi.samples

#### Summary of estimation ###

W.obs.ord.stats <- summary_stats(W.ord.post, data$true.W.obs.ord)
beta.stats <- summary_stats(beta.post, data$true.beta)
Sigma.stats <- summary_stats(Sigma.post, data$true.Sigma)
phi.stats <- summary_stats(phi.post, data$true.phi)

summary_values <- list(
  W.obs.ord.stats = W.obs.ord.stats,
  beta.stats = beta.stats,
  Sigma.stats = Sigma.stats,
  phi.stats = phi.stats
)

list.save(summary_values, file = "summary.Rdata")
#### ESS calculation for component chains ####

beta.ess <- compute_ess(beta.post)
Sigma.ess <- compute_ess(Sigma.post)
phi.ess <- compute_ess(phi.post)
W.obs.ord.ess <- compute_ess(W.ord.post)

ess_values <- list(
  beta.ess = beta.ess,
  Sigma.ess = Sigma.ess,
  phi.ess = phi.ess,
  W.obs.ord.ess = W.obs.ord.ess
)
# Saving ESS
list.save(ess_values, file = "ess.Rdata")

pred.locs <- data$pred.locs
N.pred <- data$N.pred
pred.ord <- max_min(pred.locs)
ord <- c(obs.ord, pred.ord + nrow(obs.locs))
locs.ord <- rbind(obs.locs, pred.locs)[ord, , drop = FALSE]
distlocs.ord <- rdist(locs.ord)
NNarray.all <- neighbor_matrix(locs.ord, m)
distlocs.nn <- dist.nn(locs.ord, NNarray.all)
NNarray.pred=NNarray.all[N.obs + (1:N.pred),]


#neighbor_matrix(pred.locs, m)

X.pred <- data$X.pred
Y.pred.true <- data$Y.pred.true

X.pred.ord <- X.pred[pred.ord, , drop = FALSE]
Y.pred.ord.true <- Y.pred.true[pred.ord, , drop = FALSE]

U.joint <- U.sgv(distlocs.nn, NNarray.all, phi, 
                 nu, m)


U.obs.col <- U.joint[ , 1:N.obs]
U.pred.col <- U.joint[ , (N.obs+1):(N.pred+N.obs)]
W.pred.prec <- crossprod.spam(U.pred.col)
W.pred.var <- solve.spam(as.spam(W.pred.prec))
chol.W.pred.var <- chol(W.pred.var)


thinned.post.samples <- thinning(pred.iters, niters, burnin, 
                                 W.ord.samples = W.ord.post, 
                                 beta.samples = beta.post , 
                                 Sigma.samples = Sigma.post,
                                 chol.Sigma.samples = chol.Sigma.post,
                                 phi.samples = phi.post)

W.obs.ord.samples <- thinned.post.samples$thinned.W.ord.samples
beta.samples <- thinned.post.samples$thinned.beta.samples
Sigma.samples <- thinned.post.samples$thinned.Sigma.samples
chol.Sigma.samples <- thinned.post.samples$thinned.chol.Sigma.samples
phi.samples <- thinned.post.samples$thinned.phi.samples


Y.pred.joint <- predictive.samples(W.obs.ord.samples, beta.samples, chol.Sigma.samples, 
                                   phi.samples, nu, m,
                                   X.obs.ord, X.pred.ord, U.joint, U.pred.col, U.obs.col, 
                                   chol.W.pred.var, N.obs, N.pred, pred.iters, q, p)



rmsp.joint <- RMSPE(Y.pred.samples = Y.pred.joint, Y.true = Y.pred.ord.true, pred.iters)

ess = list.load(file = "ess.Rdata")
