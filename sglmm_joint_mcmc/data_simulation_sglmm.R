rm(list = ls())
#
# mydir <- "C:/Users/Arghya/OneDrive - IIT Kanpur/arghya/IIT Kanpur PhD documents/Spatial and MCMC Research/Project1/sglmm"

library(this.path)

mydir <- this.path::here()
setwd(mydir)

#mydir_lab <- "C:/Users/HP/OneDrive - IIT Kanpur/arghya/IIT Kanpur PhD documents/Spatial and MCMC Research/Project1/Computing/Multivariate GP_ncp"
#setwd(mydir_lab)

library(fields)
library(fBasics)  # For vectorize columns of a matrix
library(rlist)

set.seed(1000)

# dimension of the random field
q <- 2


N <- 2e2 # No. of spatial locations

#simulating N locations over [0,1]^2 in locations matrix

locations <- matrix(0, nrow = N, ncol = 2)

locations[, 1] <- runif(N)
locations[, 2] <- runif(N)

# Adding a mean term

X = cbind(1, locations)

#design matrix in concatenated form of dim (Np * 1)
X.vec <- vec(t(X))

#Number of features in the spatial model
p <- ncol(X)

# True value of regression matrix beta
true.beta <- matrix(rnorm(p * q, mean = 0, sd = 2), nrow = p, ncol = q)

#beta.mat is the kronecker product of the coefficient matrix across locations

beta.mat <- diag(N) %x% t(true.beta)

# Fixed effect (I_n * B^T) X
mu.vec <- c(beta.mat %*% X.vec)


true.Sigma <- matrix(c(3, 1, 1, 2), nrow = q, ncol = q, byrow = TRUE)

# true.phi

true.phi <- 0.1

# true.nu

true.nu <- 0.5

#true.r

true.r <- 0.8

# Calculates the Euclidean distance matrix

distmat <- rdist(locations)

# Calculates the correlation matrix from matern kernel

K <- Matern(distmat, range = true.phi, smoothness = true.nu)

# Calculates the separable covariance matrix with nugget effect

Omega <- (true.r*K  + (1-true.r)*diag(N)) %x% true.Sigma

# Generating the response vector of dimension Nq*1

W.vec <- mu.vec + c(t(chol(Omega)) %*% rnorm(N*q))
W.mat <- array(W.vec, dim = c(N, q))

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



Y.vec <- rep(0, N*q)  

for(i in 1: N*q){
  if(i %% q == 1){
    Y.vec[i] = rpois(1, lambda = inv.link(W.vec[i], link = "exp"))
  }else{
    Y.vec[i] = rpois(1, lambda = inv.link(W.vec[i], link = "exp"))
  }
}

Y <- array(Y.vec, dim = c(N, q))

# Saving necessary parameters and data (all in matrix form)
save(N, p, q, locations, X, Y, W.mat, true.beta,
     true.Sigma, true.phi, true.nu, true.r, file = "sglmm.mgp.simulation.Rdata")

