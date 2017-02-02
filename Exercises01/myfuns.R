###########################################################
######### Created by Spencer Woody on 31 Jan 2017 #########
###########################################################

my.lm <- function(X, y) {
	#  Custom function for linear regression 
	#  
	#  Note: this function assumes that X already has an intercept term
	# (or doesn't, if we want to force OLS through the origin)
	#
	# INPUTS:
	# X is the design matrix
	# y is the response vector
	#
	# OUTPUTS
	# a list of...
	#  Beta.hat is a vector of estimates of the coefficients
	#  Beta.SE is a vector of the standard errors of the coefficients
	#  Beta.t is a vector of t-scores of the coefficients
	#  Beta.p is the p-value for each coefficient
	#  RSS is the residual sum of squares
	#  Var.hat is the estimated variance of homoscedastic residuals
	#  R.sq is the R-squared value
	#  R.sqadj is the adjusted R-squared value
	#  
	
	N <- nrow(X)
	p <- ncol(X)
	
	XtX <- crossprod(X)
	
	# Calculate beta.hat
	beta.hat <- solve(XtX, crossprod(X, y))
	
	# Calculate predicted values and residuals
	y.hat <- crossprod(t(X), beta.hat)
	res <- y - y.hat
	
	rss <- sum(res^2)
	
	# Calculate \hat{sigma^2}
	var.hat <- rss / (N - p)
	
	# Calculate covariance matrix of beta and SE's of beta
	var.beta <- var.hat * solve(XtX)
	beta.SE <- diag(var.beta) ^ 0.5
	
	# Calculate t-score of each beta
	beta.t <- beta.hat / beta.SE
	
	# Calculate p-values for coefficients
	beta.p <- 2 * (1 - pt(abs(beta.t), N - p))
	
	# Calculate r-squared and adjusted r-squared
	r.sq <- 1 - rss / sum((y - mean(y))^2)
	r.sqadj <- r.sq - (1 - r.sq) * (p - 1) / (N - p - 2)
	
	# Create a list of calculated values, return it back
	mylist <- list(Beta.hat = beta.hat, Beta.SE = beta.SE, 
				   Beta.t = beta.t, Beta.p = beta.p, RSS = rss, Var.hat = var.hat,
				   R.sq = r.sq, R.sqadj = r.sqadj, Res = res)
	return(mylist)
}

microbenchmark(sqrt(1e9), 1e9^0.5)

my.boot <- function(X, y, NN = 10000){
	#  Give bootstrapped estimate of covariance matrix of betas 
	#  Note: this function assumes that X already has an intercept term
	# (or doesn't, if we want to force OLS through the origin)
	#
	# INPUTS:
	# X is the design matrix
	# y is the response vector
	# N is the number of bootstrap simulations
	#
	# OUTPUTS: 
	#
	#
	#
	#
	
	N <- nrow(X)
	p <- ncol(X)
	
	XtX    <- crossprod(X)
	XtXinv <- solve(XtX)
	
	# Calculate beta.hat
	beta.hat <- solve(XtX, crossprod(X, y))
	
	# Calculate predicted values and residuals
	y.hat <- crossprod(t(X), beta.hat)
	res <- y - y.hat
	
	var.hat <- sum(res^2) / (N - p)
	
	# Run bootstrap
	var.hat.star <- c()
	beta.star <- matrix(nrow = p, ncol = NN)
	
	for(i in 1:NN) {
		sample.i <- sample(1:N, N, replace = T)
		res.star <- res[sample.i]
		var.hat.star <- c(var.hat.star, sum(res.star^2) / (N - p))
		
		y.star <- y.hat + res.star
		
		beta.star[, i] <- crossprod(XtXinv, crossprod(X, y.star))
	}
	
	mylist <- list('XtXinv' = XtXinv, 
				   'var.hat' = var.hat, 
				   'var.hat.star' = var.hat.star)
	
	return(mylist)
}

loglik <- function(X = NULL, y = NULL, params = NULL) {
	return(TRUE)
}