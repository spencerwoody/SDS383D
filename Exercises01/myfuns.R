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

my.boot1 <- function(X, y, B = 10000){
	#  Give bootstrapped estimate of covariance matrix of betas by
	#  SAMLING **RESIDUALS**
	#  Note: this function assumes that X already has an intercept term
	# (or doesn't, if we want to force OLS through the origin)
	#
	# INPUTS:
	# X is the design matrix
	# y is the response vector
	# N is the number of bootstrap simulations
	#
	# OUTPUT: 
	# cov.star is the estimated covariance matrix of beta-hat 
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
	
	# Run bootstrap
	beta.star <- matrix(nrow = B, ncol = p)
	
	for(i in 1:B) {
		sample.i <- sample(1:N, N, replace = T)
		res.star <- res[sample.i]
		
		y.star <- y.hat + res.star
		
		beta.star[i, ] <- crossprod(XtXinv, crossprod(X, y.star))
	}
	
	cov.star <- cov(beta.star)
	
	return(cov.star)
}

my.boot2 <- function(X, y, B = 10000){
	#  Give bootstrapped estimate of covariance matrix of betas by
	#  SAMLING **POINTS x & y**
	#  Note: this function assumes that X already has an intercept term
	# (or doesn't, if we want to force OLS through the origin)
	#
	# INPUTS:
	# X is the design matrix
	# y is the response vector
	# N is the number of bootstrap simulations
	#
	# OUTPUT: 
	# cov.star is the estimated covariance matrix of beta-hat 
	#
	#
	#
	
	N <- nrow(X)
	p <- ncol(X)

	# Run bootstrap
	beta.star <- matrix(nrow = B, ncol = p)
	
	for(i in 1:B) {
		sample.i <- sample(1:N, N, replace = T)
		
		X.star <- X[sample.i, ]
		y.star <- y[sample.i, ]
		
		XtX.star <- crossprod(X.star)
		
		beta.star[i, ] <- solve(XtX.star, crossprod(X.star, y.star))
	}
	
	cov.star <- cov(beta.star)
	
	return(cov.star)
}

my.mvn <- function(n, mu, Sigma) {
	#  Simulate n draws from MVN(mu, Sigma)
	#  
	#  Note: this function assumes that X already has an intercept term
	# (or doesn't, if we want to force OLS through the origin)
	#
	# INPUTS:
	# n is the number of draws
	# mu is the mean vector
	# Sigma is the covariance matrix
	#
	# OUTPUT: 
  	# x is matrix of n draws from MVN(mu, Sigma) [with n rows, p columns]
	#
	
	p <- length(mu)
	
	if ((nrow(Sigma) != ncol(Sigma)) | (det(Sigma) <= 0) | (nrow(Sigma) != p) ) {
		return("Try again...")
	}
	
	z     <- matrix(rnorm(n*p), nrow = p)
	mumat <- matrix(rep(mu, n), nrow = p)
	
	Lt <- chol(Sigma)
	
	x <- crossprod(Lt, z) + mumat
	
	return(t(x))
}

loglik <- function(X = NULL, y = NULL, params = NULL) {
	return(TRUE)
}