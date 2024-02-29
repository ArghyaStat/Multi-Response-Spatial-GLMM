
# Data Simulation

# Simulating with  Phi = 0.1

library(devtools)
library(Rcpp)
library(fields)
library(plot3D)
library(geoR)

N <- 2e2 # No. of spatial locations 

locations <- matrix(0, nrow = N, ncol = 2)

locations[, 1] <- runif(N)
locations[, 2] <- runif(N) #simulating 100 locations


plot(locations)
true_phi <- 0.1   #true value of phi


d <- rdist(locations) # Calculates the Euclidean distance matrix
cov <- exp(-d/true_phi) # Calculates the correlation matrix
Y <- c(t(chol(cov)) %*% rnorm(N)) # Draws a sample from multivariate normal


# Step 1: Save the C++ code to a file named source.cpp

# Step 2: Compile the C++ code


sourceCpp("mcmc_with_phi.cpp")

# Set the parameters
init_phi <- 0.5
nu <- 0.5
niters <- 1000
h <- 0.05




# Run the Metropolis-Hastings algorithm
phi_chain <- metropolis_hastings(init_phi = 0.5, nu = 1, niters = 1000, h = 0.05, Y = Y, locations = locations)
