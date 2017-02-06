###########################################################
######### Created by Spencer Woody on 29 Jan 2017 #########
###########################################################

library(microbenchmark)
library(ggplot2)
library(mlbench)
library(xtable)
library(mvtnorm)

source("myfuns.R")

###########################################################
#### Linear regression

# Import the data and remove missing values
ozone = data(Ozone, package='mlbench')
ozone = na.omit(Ozone)[,4:13]

# Create response vector and design matrix (with intercept)
y <- as.matrix(ozone[,1])
X <- as.matrix(ozone[,2:10])

N   <- nrow(X)
int <- rep(1, N)
X   <- cbind(int, X)

microbenchmark(
	model1 <- lm(formula = y ~ X - 1), 
	model2 <- my.lm(X, y)
	)
# my code runs about six times as fast :)

summary(model1)

model2$Beta.hat
model2$Beta.SE
model2$Beta.t 
model2$Beta.p 

###########################################################
#### Bootstrapping

# Bootstrap estimate of covariance matrix of sampling distribution 
# of betahat, resampling residuals
my.cov1 <- my.boot.res(X, y)
xtable(my.cov1, display = rep("e", 11), digits = 2)

# Bootstrap estimate of covariance matrix of sampling distribution 
# of betahat, resampling pairs x and y
my.cov2 <- my.boot.pairs(X, y)
xtable(my.cov2, display = rep("e", 11), digits = 2)

# Parametric estimate of covariance matrix of sampling distribution of betahat
cov.para <- model2$Var.hat * solve(crossprod(X))
xtable(cov.para, display = rep("e", 11), digits = 2)

# Generate MVN random variables 
en <- 50
mu <- c(3, -2)
Sigma <- matrix(c(1, 0.5, 0.5, 2), nrow = 2)

# MLE estimation of mean vector and covariance matrix 
v <- my.mvn(en, mu, Sigma)
v.mle <- mle.mvn(v)

# Bootstrap the sampling distribution of MLEs
BOOT <- my.boot.mle(v)

# Plot bootstrap results of mean vector 
pdf("mu.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), din = c(10,5), family = "serif")
hist(BOOT$mu.boot[, 1], 
	main = sprintf("Sampling distribution of mu[1], N = %i", en),
	xlab = "Boostrap estimate of mu[1]", 
	freq = T)
abline(v = mu[1], col = "red")
abline(v = v.mle$mu.hat[1], col = "blue", lty = 2)
hist(BOOT$mu.boot[, 2], 
	main = sprintf("Sampling distribution of mu[2], N = %i", en),
	xlab = "Boostrap estimate of mu[2]", 
	freq = T)
abline(v = mu[2], col = "red")
abline(v = v.mle$mu.hat[2], col = "blue", lty = 2)
dev.off()

# Plot bootstrap results of covariance matrix
pdf("Sigma.pdf", width = 10, height = 10)
par(mfrow = c(2, 2), family = "serif")
hist(BOOT$Sigma.boot[, 1], 
	main = sprintf("Sampling distribution of cov[x1, x1], N = %i", en),
	xlab = "Boostrap estimate of cov[x1, x1]", 
	freq = T)
abline(v = Sigma[1, 1], col = "red")
abline(v = v.mle$Sigma.hat[1, 1], col = "blue", lty = 2)
plot.new()
hist(BOOT$Sigma.boot[, 3], 
	main = sprintf("Sampling distribution of cov[x1, x2], N = %i", en),
	xlab = "Boostrap estimate of cov[x1, x1]", 
	freq = T)
abline(v = Sigma[2, 1], col = "red")
abline(v = v.mle$Sigma.hat[2, 1], col = "blue", lty = 2)
hist(BOOT$Sigma.boot[, 2], 
	main = sprintf("Sampling distribution of cov[x2, x2], N = %i", en),
	xlab = "Boostrap estimate of cov[x2, x2]", 
	freq = T)
abline(v = Sigma[2, 2], col = "red")
abline(v = v.mle$Sigma.hat[2, 2], col = "blue", lty = 2)
dev.off()

# Example of triangular matrices
A <- matrix(c(1,2,3,4,5,6,7,8,9), nrow = 3)
A[lower.tri(A)]