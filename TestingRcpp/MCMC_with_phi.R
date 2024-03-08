
# Data Simulation

# Simulating with  Phi = 0.1

library(devtools)
library(Rcpp)
library(fields)
library(plot3D)
library(geoR)

N <- 100 # No. of spatial locations 

locations <- matrix(0, nrow = N, ncol = 2)

locations[, 1] <- runif(N)
locations[, 2] <- runif(N) #simulating N locations


plot(locations)
true_phi <- 0.9   #true value of phi


d <- rdist(locations) # Calculates the Euclidean distance matrix
cov <- exp(-d/true_phi) # Calculates the correlation matrix
Y <- c(t(chol(cov)) %*% rnorm(N)) # Draws a sample from multivariate normal


# Calling the C++ code from a file named mcmc_with_phi.cpp


sourceCpp("mcmc_with_phi_new.cpp")

# Set the parameters
init_phi <- 0.5
nu <- 0.5
niters <- 2
h <- 0.3

# tem <- function(locations, phi)
# {
#   n.l <- dim(locations)[1]
#   out.mat <- matrix(0, nrow = n.l, ncol = n.l)
#   for(i in 1:n.l)
#   {
#     for(j in 1:n.l)
#     {
#       xi <- locations[i, ]
#       xj <- locations[j, ]
#       h <-  sqrt(sum(( xi - xj)^2))
#       out.mat[i,j] <- exp(-h/phi)
#     }
#   }
#   return(out.mat)
# }

# matern(dist = .2, phi = 0.5, nu = .5)
# foo <- cov_mat(locations = locations,  phi = .5, nu = .5)
# foo - tem(locations, phi = .5)
# eigen(foo, only.values = TRUE)$values
# 
# nu <- .5
# phi <- .5
# h <- 2
# temp <- ((sqrt(2*nu)*h)/phi)^nu * 2^(1 - nu)/gamma(nu)
# exp(-h/phi)/temp
# 
# matern(dist = h, phi = phi, nu = nu)


#bessel_k(dist, nu, 1)
# Run the Metropolis-Hastings algorithm
phi_chain <- metropolis_hastings(init_phi = 0.5, nu = 0.5, niters = 1e3, h = 0.05, Y = Y, locations = locations)
plot.ts(phi_chain)
