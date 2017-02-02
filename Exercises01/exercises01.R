###########################################################
######### Created by Spencer Woody on 29 Jan 2017 #########
###########################################################

library(microbenchmark)
library(ggplot2)
library(mlbench)

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

BOOT <- my.boot(X, y)