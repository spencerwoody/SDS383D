###########################################################
######### Created by Spencer Woody on 29 Jan 2017 #########
###########################################################

library(microbenchmark)
library(ggplot2)
library(mlbench)
library(xtable)

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
my.cov1 <- my.boot1(X, y)
xtable(my.cov1, display = rep("e", 11), digits = 2)

# Bootstrap estimate of covariance matrix of sampling distribution 
# of betahat, resampling pairs x and y
my.cov2 <- my.boot2(X, y)
xtable(my.cov2, display = rep("e", 11), digits = 2)

# Parametric estimate of covariance matrix of sampling distribution of betahat
cov.para <- model2$Var.hat * solve(crossprod(X))
xtable(cov.para, display = rep("e", 11), digits = 2)


# MLE estimation of mean vector and covariance matrix 