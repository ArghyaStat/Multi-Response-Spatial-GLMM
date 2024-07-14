
### Predictive Analysis




W.tilde.predictive <- function(W.tilde.obs, X.obs, X.pred, beta, Sigma, phi, nu, r, distmat.joint){
  
  N.obs <- nrow(X.obs)
  N.pred <- nrow(X.pred)
  
  K.tilde.joint <- r* Matern(distmat.joint, range = phi, smoothness = nu) + (1-r)* diag(N.pred+N.obs)
  
  K.tilde.pred <- K.tilde.joint[1:N.pred, 1:N.pred]
  K.tilde.obs <-  K.tilde.joint[(N.pred+1): (N.pred+N.obs), (N.pred+1):(N.pred+N.obs)]
  K.tilde.cross <- K.tilde.joint[(1:N.pred), (N.pred+1):(N.pred+N.obs)]
  
  chol.K.tilde.obs <- chol(K.tilde.obs)
  K.tilde.obs.inv <- chol2inv(chol.K.tilde.obs)
  K.tilde.reg <- K.tilde.cross %*% K.tilde.obs.inv
  
  
  
  W.tilde.pred.mean <- X.pred %*% beta + K.tilde.reg %*% (W.tilde.obs - X.obs %*% beta)
  
  K.tilde.cond <- K.tilde.pred - K.tilde.reg %*% t(K.tilde.cross)
  chol.K.tilde.cond <- chol(K.tilde.cond)
  
  chol.Sigma <- chol(Sigma)
  
  W.pred.samp <- matrix(rnorm(N.pred*q, 0, 1), nrow = N.pred , ncol = q)
  
  W.tilde.pred <- W.tilde.pred.mean + t(chol.K.tilde.cond) %*% W.pred.samp %*% chol.Sigma
  
  W.tilde.pred
  
}

Y.predictive <- function(W.tilde.pred){
  
  N.pred <- nrow(W.tilde.pred)
  q <- ncol(W.tilde.pred)
  
  Y.mat <- matrix(NA, nrow = N.pred, ncol = q)
  
  
  for(i in 1:N.pred){

    #Y.mat[i,1] <- rnorm(1, mean = W.tilde.pred[i,1], sd = 1)
    Y.mat[i,1] <- rpois(1, lambda = exp(W.tilde.pred[i,1]))

  }
  
  Y.mat
}
  

thinning <- function(pred.iters, niters, burnin,
                     W.tilde.samples, beta.samples, Sigma.samples, phi.samples, nu.samples, r.samples){
  
  
  
  if ((niters-burnin) %% pred.iters != 0) {
    stop("niters must be divisible by pred.iters")
  }
  
  thin.block <- (niters-burnin)/pred.iters
  thinned.indices <- seq(burnin + 1, niters, by = thin.block)
  
  thinned.W.tilde.samples <- list()
  thinned.beta.samples <- list()
  thinned.Sigma.samples <- list()
  
  
  for (i in thinned.indices) {
    
    thinned.W.tilde.samples <- c(thinned.W.tilde.samples, list(W.tilde.samples[[i]]))
    thinned.beta.samples <- c(thinned.beta.samples, list(beta.samples[[i]]))
    thinned.Sigma.samples <- c(thinned.Sigma.samples, list(Sigma.samples[[i]]))
    
  }
  
  thinned.phi.samples <- phi.samples[thinned.indices]
  thinned.nu.samples <- nu.samples[thinned.indices]
  thinned.r.samples <- r.samples[thinned.indices]
  
  thinned.samples <- list("thinned.W.tilde.samples" = thinned.W.tilde.samples,
                        "thinned.beta.thinned.samples" = thinned.beta.samples,
                        "thinned.Sigma.samples" = thinned.Sigma.samples,
                        "thinned.phi.samples" = thinned.phi.samples,
                        "thinned.nu.samples" = thinned.nu.samples,
                        "thinned.r.samples" = thinned.r.samples)
  
  return(thinned.samples)
  
}


predictive.samples <- function(pred.iters, X.obs, X.pred, distmat.joint,
                         W.tilde.samples, beta.samples, Sigma.samples, phi.samples, nu.samples, r.samples){
  
  Y.pred <- replicate(pred.iters, matrix(0, N.pred, q), simplify = F)
  W.tilde.pred <- replicate(pred.iters, matrix(0, N.pred, q), simplify = F)
  
  for(i in 1:pred.iters){
    
    W.tilde.pred[[i]] <-  W.tilde.predictive(W.tilde.samples[[i]], X.obs, X.pred, 
                                             beta.samples[[i]], 
                                             Sigma.samples[[i]], 
                                             phi.samples[[i]], 
                                             nu.samples[[i]], 
                                             r.samples[[i]], distmat.joint)
    Y.pred[[i]] <- Y.predictive(W.tilde.pred[[i]])
  }
  
  Y.pred
  
}


load("wildfire.Rdata")
theta.chain <- list.load(file = "wildfire.MCMC.chain.Rdata")

distmat.joint <- rdist(joint.loc)

niters <- 1e5
burnin <- 2e4
pred.iters <- 8e4

W.tilde.post <- theta.chain$W.tilde.sample
beta.post <- theta.chain$beta.sample
Sigma.post <- theta.chain$Sigma.sample
phi.post <- theta.chain$phi.sample
nu.post <- theta.chain$nu.sample
r.post <- theta.chain$r.sample

thinned.post.samples <- thinning(pred.iters, niters, burnin, 
                                 W.tilde.samples = W.tilde.post, 
                                 beta.samples = beta.post , 
                                 Sigma.samples = Sigma.post,
                                 phi.samples = phi.post, 
                                 nu.samples = nu.post, 
                                 r.samples = r.post)

W.tilde.samples <- thinned.post.samples$thinned.W.tilde.samples
beta.samples <- thinned.post.samples$thinned.beta.thinned.samples
Sigma.samples <- thinned.post.samples$thinned.Sigma.samples
phi.samples <- thinned.post.samples$thinned.phi.samples
nu.samples <- thinned.post.samples$thinned.nu.samples
r.samples <- thinned.post.samples$thinned.r.samples


par(mfrow = c(1,1))

pdf("traceplot_phi.pdf", width=5 ,height=5)
trace.phi <- plot.ts(phi.samples, ylab = "phi", main = "Traceplot of phi")
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
    plot(1:pred.iters, sapply(beta.samples, function(x) x[i, j]), 
         type='l', col=1, xlab='Iteration', ylab='beta',
         main=labels.beta[(i-1)*q + j])
    
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
    plot(1:pred.iters, sapply(Sigma.samples, function(x) x[i, j]), 
         type='l', col=1, xlab='Iteration', ylab='Sigma',
         main=labels.Sigma[(i-1)*q + j])
    
  }
}
dev.off()

extremes.W.tilde <- function(W.tilde.samples){
  
  thinned.N <- length(W.tilde.samples)
  
  out <- matrix(NA, nrow = thinned.N, ncol = 2)
  
  for(i in 1: thinned.N){
    
    out[i,1] <- min(W.tilde.samples[[i]][,2])
    out[i,2] <- max(W.tilde.samples[[i]][,2])
  }
  
  return(out)
}



extremes.W.tilde.count <- extremes.W.tilde(W.tilde.samples = W.tilde.samples)

save(extremes.W.tilde.count, file = "extremes.W.tilde.count.Rdata")

Y.pred.joint <- predictive.samples(pred.iters, X.obs, X.pred, distmat.joint,
                                   W.tilde.samples, beta.samples, Sigma.samples, phi.samples, nu.samples, r.samples)

Y.pred.sep <- predictive.samples(pred.iters, X.obs, X.pred, distmat.joint,
                                      W.tilde.samples, beta.samples, Sigma.samples, phi.samples, nu.samples, r.samples)

Y.pred.joint.mean <- Reduce("+", Y.pred.joint)/pred.iters
Y.pred.sep.mean <- Reduce("+", Y.pred.sep)/pred.iters

RMSP.predict <- function(Y.pred.samples, Y.true, pred.iters){
  
  out <- 0
  
  q <- ncol(Y.true)
  
  for(i in 1:pred.iters){
    
    for(j in 1:q){
      
      out <- out + sum((Y.pred.samples[[i]][,q] - Y.true[,q])^2)
    }
    
  }
  
  return(sqrt(out/pred.iters))
}

rmsp.joint <- RMSP.predict(Y.pred.samples = Y.pred.joint, Y.true = Y.pred.true, pred.iters = 8e4)

mse.joint <- norm(Y.pred.joint.mean, Y.pred.true, type = "2")
mse.sep <- norm(Y.pred.sep.mean, Y.pred.true, type = "2")


# Spatial Visualization

# Load necessary libraries
library(ggplot2)
library(viridis)
library(gridExtra)

# Function to plot the two types of predicted responses along with their true values
plot_responses <- function(Y.pred.true, Y.pred.est, pred.loc) {
  # Create data frames for true and predicted values
  true_data <- data.frame(x = pred.loc[,1], y = pred.loc[,2], continuous = Y.pred.true[,1], count = Y.pred.true[,2])
  pred_data <- data.frame(x = pred.loc[,1], y = pred.loc[,2], continuous = Y.pred.est[,1], count = Y.pred.est[,2])
  
  
  scale.cont <- range(c(true_data$continuous, pred_data$continous))
  scale.count <- range(c(true_data$count, pred_data$count))
  # True continuous response plot
  p1 <- ggplot(true_data, aes(x = x, y = y)) +
    geom_point(aes(color = continuous), size = 2) +
    scale_color_viridis_c(name = "Continuous", limits = scale.cont) +
    labs(title = "True Continuous",
         x = "Latitute",
         y = "Longitude") +
    theme_minimal()
  
  # Predicted continuous response plot
  p2 <- ggplot(pred_data, aes(x = x, y = y)) +
    geom_point(aes(color = continuous), size = 2) +
    scale_color_viridis_c(name = "Continuous", limits = scale.cont) +
    labs(title = "Predicted Continuous",
         x = "Latitute",
         y = "Longitude") +
    theme_minimal()
  
  # True count response plot
  p3 <- ggplot(true_data, aes(x = x, y = y)) +
    geom_point(aes(color = count), size = 2) +
    scale_color_viridis_c(name = "Count", limits = scale.count) +
    labs(title = "True Count",
         x = "Latitute",
         y = "Longitude") +
    theme_minimal()
  
  # Predicted count response plot
  p4 <- ggplot(pred_data, aes(x = x, y = y)) +
    geom_point(aes(color = count), size = 2) +
    scale_color_viridis_c(name = "Count", limits = scale.count) +
    labs(title = "Predicted Count",
         x = "Latitute",
         y = "Longitude") +
    theme_minimal()
  
  # Arrange the plots in a 2x2 grid
  grid.arrange(p3, p4, nrow = 2, top = "True and Predicted Responses")
}



# Call the function with the example data

pdf("comparative_visualization_joint.pdf", width = 7.5 ,height = 6)
joint.predict.plot <- plot_responses(Y.pred.true = Y.pred.true, Y.pred.est = Y.pred.joint.mean, pred.loc = pred.loc)
dev.off()

pdf("comparative_visualization_marginal.pdf", width = 7.5 ,height = 6)
sep.predict.plot <- plot_responses(Y.pred.true = Y.pred.true, Y.pred.est = Y.pred.sep.mean, pred.loc = pred.loc)
dev.off()


