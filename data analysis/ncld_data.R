rm(list = ls())


library(this.path)

mydir <- this.path::here()
setwd(mydir)

library(fields)
library(ggplot2)
library(viridis)
library(plot3D)
library(fBasics)  # For vectorize columns of a matrix
library(MCMCpack)
library(mvtnorm)

load(file.path("cleaned_data_california.Rdata"), data <- new.env())
data <- as.list(data)

lat <- data$latS
lon <- data$lonS

locations <- cbind(lat, lon)

# Adding a mean term

X <- cbind(1, locations)

#Number of features in the spatial model
p <- ncol(X)

# Response 

Y_BA <- log(1+data$BA.mat[,19])
Y_CNT <- data$CNT.mat[,19]
  
#log(1+data$CNT.mat[,19])
  
  

  
#log(1+data$CNT.mat[,19])

# q <- 1
# 
Y <- as.matrix(Y_CNT)

#Y <- cbind(Y_BA, Y_CNT)

N <- nrow(Y)
q <- ncol(Y)

#### Data partitioning

N.pred <- ceiling(0.2*N)
N.obs <- N - N.pred

pred.indices <- sample(1:N, N.pred)
pred.loc <- locations[pred.indices, ]
obs.loc <- locations[-pred.indices, ]

distmat.obs <- rdist(obs.loc)
diameter <- max(distmat.obs)

joint.loc <- rbind(pred.loc, obs.loc)

X.obs <- X[-pred.indices,]
X.pred <- X[pred.indices,]

Y.pred.true <- as.matrix(Y[pred.indices,])
Y.obs <- as.matrix(Y[-pred.indices,])



# Saving necessary parameters and data (all in matrix form)
save(N, N.obs, N.pred, p, q, joint.loc, obs.loc, pred.loc, X.obs, X.pred, 
     diameter, Y.obs, Y.pred.true, file = "wildfire.Rdata")



### Exploratory Data Analysis

lm.BA <- lm(Y_BA ~ locations)
lm.CNT <- lm(log(1+Y_CNT) ~ locations)


resid.BA <- lm.BA$residuals
resid.CNT <- lm.CNT$residuals

# Histograms of residuals of responses removing the effect of covariates



pdf("Hist_resid_log(1+BA)_and_log(1+CNT).pdf", width=10 ,height=5)


hist1 <- hist(resid.BA, breaks=10, probability=TRUE, 
              main="Histogram of residuals of log(1+BA)", xlim = c(-6,6),
              xlab="Value", ylab="Empirical density")
dev.off()

pdf("Hist_resid_log(1+CNT).pdf", width = 5, height = 5)

hist2 <- hist(resid.CNT, breaks=10, probability=TRUE, 
              main="Histogram of residuals of log(1+CNT)",
              xlab="Value", ylab="Empirical density")
dev.off()





library(sp)
library(gstat)

#### Variogram analysis of BA

df.BA <- data.frame(lat = locations[ , 1], lon = locations[ , 2], values = resid.BA)
coordinates(df.BA) <- ~ lat+lon
emp.varigram.BA <- variogram(values~1, df.BA)

fitted.varigram.BA <- fit.variogram(emp.varigram.BA, vgm(psill = var(emp.varigram.BA$gamma), "Mat", nugget = NA, fit.kappa = TRUE))


fit.BA <- vgm(psill = fitted.varigram.BA$psill[2], "Mat", range = fitted.varigram.BA$range[2], kappa = 0.5)

est.sigmasq.BA <- fitted.varigram.BA$psill[2]
est.phi.BA <- fitted.varigram.BA$range[2]

#### Variogram analysis of CNT


df.CNT <- data.frame(lat = locations[ , 1], lon = locations[ , 2], values = resid.CNT)
coordinates(df.CNT) <- ~ lat+lon
emp.varigram.CNT <- variogram(values~1, df.CNT)



fitted.varigram.CNT <- fit.variogram(emp.varigram.CNT, vgm(psill = var(emp.varigram.CNT$gamma), "Mat", nugget = NA, fit.kappa = TRUE))


fit.CNT <- vgm(psill = fitted.varigram.CNT$psill[2], "Mat", range = fitted.varigram.CNT$range[2], kappa = 0.5)

est.sigmasq.CNT <- fitted.varigram.CNT$psill[2]
est.phi.CNT <- fitted.varigram.CNT$range[2]


# Variograms of residuals of responses removing the effect of covariates

pdf("Variogram_resid_log(1+BA).pdf", width=5 ,height=5)

plot(emp.varigram.BA, fit.BA, pch = 19, lwd = 2, col = "black", main = "Variogram of residuals of log(1+BA)")

dev.off()

pdf("Variogram_resid_log(1+CNT).pdf", width = 5, height = 5)

plot(emp.varigram.CNT, fit.CNT, pch = 19, lwd = 2, col = "red", main = "Variogram of residuals of log(1+CNT)")

dev.off()
