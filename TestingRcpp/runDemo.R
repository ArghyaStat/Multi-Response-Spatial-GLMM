##################################
## Main interactive code
##################################
library(Rcpp)
sourceCpp("demo.cpp")
# code for AR(1)
ar1 <- function(n = 1e3, rho = .5, start = 0)
{
  output <- ar1Cpp(n = n, rho = rho, start = start)
  return(output)
}

slowar1 <- function(n = 1e3, rho = .5, start = 0)
{
  out <- numeric(length = n)
  out[1] <- start
  errors <- rnorm(n)
  for(i in 2:n)
  {
    out[i] <- rho*out[i-1] + errors[i]
  }
  return(out)
}

chain <- slowar1(rho = .99, star = 150)
plot.ts(chain)

library(rbenchmark)
benchmark(ar1(n = 1e6), slowar1(n = 1e6), replications = 10)


