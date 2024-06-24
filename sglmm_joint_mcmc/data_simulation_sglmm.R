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

set.seed(1e3)

# dimension of the random field
q <- 2

N <- 2e2 # No. of spatial locations

#simulating N locations over [0,1]^2 in locations matrix

N.pred <- ceiling(0.2*N)
N.obs <- N - N.pred

locations <- matrix(0, nrow = N, ncol = 2)

locations[, 1] <- runif(N)
locations[, 2] <- runif(N)

pred.indices <- sample(1:N, N.pred)
pred.loc <- locations[pred.indices, ]
obs.loc <- locations[-pred.indices, ]

joint.loc <- rbind(pred.loc, obs.loc)

# Adding a mean term

X.obs <- cbind(1, obs.loc)
X.pred <- cbind(1, pred.loc)

X <- cbind(1, locations)

#Number of features in the spatial model
p <- ncol(X.obs)

# True value of regression matrix beta
true.beta <- matrix(rnorm(p * q, mean = 0, sd = 2), nrow = p, ncol = q)


# Fixed effect 
mu.vec <- vec(X %*% true.beta)


true.Sigma <- matrix(c(3, 2, 2, 4), nrow = q, ncol = q, byrow = TRUE)

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

Omega <- true.Sigma %x% (true.r*K  + (1-true.r)*diag(N))

# Generating the response vector of dimension Nq*1

true.W.tilde.vec <- mu.vec + c(t(chol(Omega)) %*% rnorm(N*q))
true.W.tilde <- array(true.W.tilde.vec, dim = c(N, q))

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

for(i in 1:N*q){
  if(i %% q == 1){
    Y.vec[i] = true.W.tilde.vec[i]
      
  }else{
    Y.vec[i] = rpois(1, lambda = inv.link(true.W.tilde.vec[i], link = "exp"))
  }
}

Y <- array(Y.vec, dim = c(N, q))
Y.pred.true <- Y[pred.indices,]
Y.obs <- Y[-pred.indices,]

true.W.tilde.obs <- true.W.tilde[-pred.indices,]

# Saving necessary parameters and data (all in matrix form)
save(N, N.obs, N.pred, p, q, joint.loc, obs.loc, pred.loc, X.obs, X.pred, 
     Y.obs, Y.pred.true, true.W.tilde.obs, true.beta, true.Sigma, true.phi, true.nu, true.r, 
     file = "sglmm.simulation.Rdata")

