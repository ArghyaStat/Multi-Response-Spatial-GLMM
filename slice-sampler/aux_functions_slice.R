
# Define function for matrix-normal density
dmatnorm.sgv <- function(X, M, chol.U.prec, chol.V.prec, NNarray, log = FALSE) {
  
  p <- nrow(X)
  q <- ncol(X)
  
  denom <- -0.5 * (p * q * log(2 * pi)) + p * sum(log(diag(chol.V.prec))) + q * sum(log(diag(chol.U.prec)))
  
  temp <- chol.U.prec %*% (X - M) 
    
  cross_mean <- tcrossprod(temp, chol.V.prec)
  
  exponent <- - 0.5 * sum(diag(crossprod(cross_mean)))
  
  if (log) {
    return(denom + exponent)
  } else {
    return(exp(denom + exponent))
  }
}


#### auto-tuning function for MH algorithm in warm up chain

tuning.update <- function(acc, att, tuning, nattempts = 50, lower = 0.8, higher = 1.2) {
  acc.rate     <- acc / att
  these.update <- att > nattempts
  these.low    <- (acc.rate < 0.2) & these.update
  these.high   <- (acc.rate > 0.3) & these.update
  
  tuning[these.low]  <- tuning[these.low] * lower
  tuning[these.high] <- tuning[these.high] * higher
  
  acc[these.update] <- 0
  att[these.update] <- 0
  
  results <- list(acc = acc, att = att, tuning = tuning)
  return(results)
}


######### Posterior Computations #########


# Conditional Log-likelihood function multi-type responsed GLM for W | Y 

log.likelihood <- function(W.obs.ord, Y.obs.ord, family) {
  # Validate inputs
  if (ncol(W.obs.ord) != ncol(Y.obs.ord)) {
    stop("The dimensions of W.obs.ord and Y.obs.ord must match.")
  }
  if (ncol(Y.obs.ord) != length(family)) {
    stop("The number of columns in Y.obs.ord must match the length of the 'family' vector.")
  }
  
  N.obs <- nrow(W.obs.ord)
  q <- ncol(W.obs.ord)
  
  like.total <- 0
  
  # Loop over observations and response types
  for (j in 1:q) {
    if (family[j] == "Gaussian") {
      # Gaussian log-likelihood: identity link
      like.total <- like.total + sum(dnorm(Y.obs.ord[, j], mean = W.obs.ord[, j], sd = 1, log = TRUE))
    } else if (family[j] == "Poisson") {
      # Poisson log-likelihood: log link (inverse is exp)
      lambda <- exp(W.obs.ord[, j])
      like.total <- like.total + sum(dpois(Y.obs.ord[, j], lambda = lambda, log = TRUE))
    } else if (family[j] == "Binomial") {
      # Binomial log-likelihood: logit link (inverse is plogis)
      prob <- plogis(W.obs.ord[, j])
      like.total <- like.total + sum(dbinom(Y.obs.ord[, j], size = 1, prob = prob, log = TRUE))
    } else if (family[j] == "Gamma") {
      # Gamma log-likelihood: log link (inverse is exp)
      shape <- exp(W.obs.ord[, j])
      rate <- 1 # Assuming fixed rate; adjust as needed
      like.total <- like.total + sum(dgamma(Y.obs.ord[, j], shape = shape, rate = rate, log = TRUE))
    } else if (family[j] == "Negative-Binomial") {
      # Negative-Binomial log-likelihood: log link (inverse is exp)
      mu <- exp(W.obs.ord[, j])
      theta <- 1 # Assuming fixed dispersion parameter; adjust as needed
      like.total <- like.total + sum(MASS::dnbinom(Y.obs.ord[, j], mu = mu, size = theta, log = TRUE))
    }else if (family[j] == "Inverse-Gaussian") {
      mu <- exp(W.obs.ord[, j])
      lambda <- 1 # Assuming fixed shape parameter; adjust as needed
      like.total <- like.total + sum(dinvgauss(Y.obs.ord[, j], mean = mu, shape = lambda, log = TRUE))
    } else if (family[j] == "Tweedie") {
      require(tweedie)
      p <- 1.5  # Assuming power parameter; adjust as needed
      phi <- 1  # Dispersion parameter; adjust as needed
      like.total <- like.total + sum(dtweedie(Y.obs.ord[, j], power = p, mu = exp(W.obs.ord[, j]), phi = phi, log = TRUE))
    } else if (family[j] == "Multinomial") {
      require(nnet)
      prob <- exp(W.obs.ord[, j]) / sum(exp(W.obs.ord[, j]))
      like.total <- like.total + sum(dmultinom(Y.obs.ord[, j], prob = prob, log = TRUE))
    } else {
      stop(paste("Unsupported family:", family[j]))
    }
    
  }
  
  return(like.total)
}


compute_ess <- function(chain) {
  # Check if the chain is a list
  if (is.list(chain)) {
    # If the chain is a list of matrices, compute ESS element-wise and keep matrix structure
    if (all(sapply(chain, is.matrix))) {
      n_iter <- length(chain)
      n_row <- nrow(chain[[1]])
      n_col <- ncol(chain[[1]])
      
      # Initialize a matrix to store ESS values for each element of the matrix chain
      ess_matrix <- matrix(NA, nrow = n_row, ncol = n_col)
      
      # Compute ESS for each element of the matrix
      for (i in 1:n_row) {
        for (j in 1:n_col) {
          # Extract the (i, j) element across all iterations and compute ESS
          element_chain <- sapply(chain, function(x) x[i, j])
          ess_matrix[i, j] <- ess(element_chain)
        }
      }
      
      return(ess_matrix)  # Return the ESS matrix
      
    } else if (all(sapply(chain, is.vector))) {
      # If the chain is a list of numeric vectors (e.g., phi.chain), compute ESS for each vector
      ess_vector <- sapply(chain, ess)  # Apply ESS to each vector in the list
      return(ess_vector)
      
    } else {
      stop("The chain contains elements that are neither matrices nor vectors.")
    }
    
  } else if (is.vector(chain)) {
    # If the chain is a single numeric vector (e.g., phi.chain), compute ESS for the entire vector
    return(ess(chain))  # Directly apply ESS to the vector
    
  } else {
    stop("Unexpected structure for chain")
  }
}




# Function to compute summary statistics including RMSE
summary_stats <- function(chain, true_value) {
  n_iter <- length(chain)
  
  if (is.matrix(chain[[1]])) {
    post.mean <- Reduce("+", chain) / n_iter
    post.median <- apply(simplify2array(chain), c(1,2), median)
    post.sd <- sqrt(Reduce("+", lapply(chain, function(x) (x - post.mean)^2)) / n_iter)
    ci.lower <- apply(simplify2array(chain), c(1,2), quantile, probs = 0.025)
    ci.upper <- apply(simplify2array(chain), c(1,2), quantile, probs = 0.975)
    coverage <- (true_value >= ci.lower) & (true_value <= ci.upper)
    rmse <- sqrt(Reduce("+", lapply(chain, function(x) (x - true_value)^2)) / n_iter)
  } else {
    post.mean <- mean(unlist(chain))
    post.median <- median(unlist(chain))
    post.sd <- sd(unlist(chain))
    ci.lower <- quantile(unlist(chain), probs = 0.025)
    ci.upper <- quantile(unlist(chain), probs = 0.975)
    coverage <- (true_value >= ci.lower) & (true_value <= ci.upper)
    rmse <- sqrt(mean((unlist(chain) - true_value)^2))
  }
  
  list(post.mean = post.mean, post.median = post.median, post.sd = post.sd, ci.lower = ci.lower,
       ci.upper = ci.upper, coverage = coverage, rmse = rmse)
}


score_function <- function(W.obs.ord, Y.obs.ord, family) {
  q <- ncol(W.obs.ord)
  n <- nrow(W.obs.ord)
  scores <- matrix(NA, nrow = n, ncol = q)
  
  for (j in 1:q) {
    W_j <- W.obs.ord[, j]
    Y_j <- Y.obs.ord[, j]
    
    if (family[j] == "Gaussian") {
      scores[, j] <- Y_j - W_j
    } else if (family[j] == "Poisson") {
      scores[, j] <- Y_j - exp(W_j)
    } else if (family[j] == "Binomial") {
      p <- exp(W_j) / (1 + exp(W_j))
      scores[, j] <- Y_j - p
    } else if (family[j] == "Gamma") {
      alpha <- 1  # Can be customized or passed
      scores[, j] <- -Y_j * exp(W_j) + alpha
    } else if (family[j] == "Negative-Binomial") {
      mu <- exp(W_j)
      theta <- 1  # Fixed dispersion
      scores[, j] <- (Y_j - mu) / (mu + mu^2 / theta)
    } else if (family[j] == "Inverse-Gaussian") {
      mu <- exp(W_j)
      lambda <- 1
      scores[, j] <- (Y_j - mu) / (mu^3 / lambda)
    } else if (family[j] == "Tweedie") {
      mu <- exp(W_j)
      p <- 1.5
      scores[, j] <- Y_j - mu^p  # Approximate; subject to link function
    } else if (family[j] == "Multinomial") {
      p <- exp(W_j) / (1 + exp(W_j))
      scores[, j] <- Y_j - p
    } else {
      stop(paste("Unsupported family:", family[j]))
    }
  }
  
  return(vec(scores))
}


