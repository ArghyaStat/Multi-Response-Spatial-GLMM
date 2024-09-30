rm(list = ls())


library(this.path)

mydir <- this.path::here()
setwd(mydir)

library(spam)
library(fields)
library(fBasics)  # For vectorize columns of a matrix
library(rlist)

set.seed(3019)

# dimension of the random field
q <- 2

N <- 2e2 # No. of spatial locations

#simulating N locations over [0,1]^2 in locations matrix


locations <- matrix(0, nrow = N, ncol = 2)

locations[, 1] <- runif(N)
locations[, 2] <- runif(N)


# Adding a mean term

X <- cbind(1, locations)

#Number of features in the spatial model
p <- ncol(X)

# True value of regression matrix beta
true.beta <- matrix(rnorm(p * q, mean = 0, sd = 1.25), nrow = p, ncol = q)


# Fixed effect 
mu.true <- X %*% true.beta


true.Sigma <- matrix(c(3,2,2,4), nrow = q, ncol = q, byrow = TRUE)

# true.phi

true.phi <- 0.1

# true.nu

true.nu <- 0.5

#true.r

true.r <- 1

# Calculates the Euclidean distance matrix

distmat <- rdist(locations)

# Calculates the correlation matrix from matern kernel

K.true <- Matern(distmat, range = true.phi, smoothness = true.nu)
K.tilde.true <- true.r*K.true + (1-true.r)*diag(N)

# Calculates the separable covariance matrix with nugget effect

chol.true.Sigma <- chol(true.Sigma)
chol.K.tilde.true <- chol(K.tilde.true)

#Omega <- true.Sigma %x% K.tilde.true

# Generating the response vector of dimension Nq*1

W.tilde.error <- matrix(rnorm(N*q, 0, 1), nrow = N , ncol = q)

true.W.tilde <- mu.true + t(chol.K.tilde.true) %*% W.tilde.error %*% chol.true.Sigma


# Some known inverse link functions

inv.link <- function(x, link) {
  if (link == "expit") {
    return(1 / (1 + exp(-x)))
  } else if (link == "exp") {
    return(exp(x))
  } else if (link == "identity") {
    return(x)
  } else if (link == "inverse") {
    return(1 / x)
  } else {
    stop("Unknown link function")
  }
}



Y <- matrix(NA, N, q)

Y[, 1] = rnorm(N, mean = true.W.tilde[, 1], sd = 1)
Y[, 2] = rpois(N, exp(true.W.tilde[, 2]))
  


N.pred <- ceiling(0.2*N)
N.obs <- N - N.pred

pred.indices <- sample(1:N, N.pred)
pred.loc <- locations[pred.indices, ]
obs.loc <- locations[-pred.indices, ]

joint.loc <- rbind(pred.loc, obs.loc)

distmat.obs <- rdist(obs.loc)
diameter <- max(distmat.obs)

X.obs <- X[-pred.indices,]
X.pred <- X[pred.indices,]

Y.pred.true <- Y[pred.indices,]
Y.obs <- Y[-pred.indices,]

true.W.tilde.obs <- true.W.tilde[-pred.indices,]

# Saving necessary parameters and data (all in matrix form)
save(N, N.obs, N.pred, p, q, joint.loc, obs.loc, pred.loc, X.obs, X.pred, diameter,
     Y.obs, Y.pred.true, true.W.tilde.obs, true.beta, true.Sigma, true.phi, true.nu, true.r, 
     file = "sglmm.bg.simulation.Rdata")

