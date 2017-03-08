###########################################################
######### Created by Spencer Woody on 11 Feb 2017 #########
###########################################################

# ===========================================================================
# Linear smoothing ==========================================================
# ===========================================================================

lin.smooth <- function(x.new, x, y, kern.fun, h) {
	# --------------------------------------------------------------------
	# Linear smoother for some kernel function 
	# --------------------------------------------------------------------
	# INPUTS:
	# x.new - a single new point for which to estimate f(x.new)
	# x - a vector of covariates from previous observations
	# y - a vector of responses from previous observations
	# kern.fun - some kernel function (e.g. Gaussian) 
	#            *** takes 2 arguments: distance (dist) and bandwidth (h)
	# h is the bandwidth for the kernel function
	# --------------------------------------------------------------------
	# OUTPUT: 
  	# weights - a vector of weights for a new observation 
	# --------------------------------------------------------------------
	
	weights <- kern.fun(dist = x - x.new, h = h) / h
	weights <- weights / sum(weights)
	
	fit <- crossprod(weights, y)
	
	return(fit)
}

kern.unif <- function(dist, h) {
	# --------------------------------------------------------------------
	# Uniform kernel function
	# --------------------------------------------------------------------
	# INPUTS:
	# dist 
	# h is
	# Sigma is the covariance matrix
	# --------------------------------------------------------------------
	# OUTPUT: 
  	# kern - the value of the uniform kernel function 
	# --------------------------------------------------------------------
	
	kern <- ( abs(dist / h) <= 1) / 2
	
	return(kern)
}

kern.norm <- function(dist, h) {
	# --------------------------------------------------------------------
	# Gaussian (normal) kernel function
	# --------------------------------------------------------------------
	# INPUTS:
	# dist the distance
	# h is the bandwidth
	# Sigma is the covariance matrix
	# --------------------------------------------------------------------
	# OUTPUT: 
  	# kern is the value of the Gaussian kernel function 
	# --------------------------------------------------------------------
	
	kern <- 1 / sqrt(2 * pi) * exp(-(dist / h)^2 / 2)
	
	return(kern)
}

make.noise <- function(x, f, res.dist, sd = NA, scale = NA, df = NA) {
	# --------------------------------------------------------------------
	# Simulate noisy response from some non-linear function
	# --------------------------------------------------------------------
	# INPUTS:
	# x - the predictor values
	# f - a function for the expected value, E(y | x) = f(x)
	# res.dist - a string for distribution of errors, either
	#            "normal" for normal errors, or 
	#            "cauchy" for cauchy error
	# sd - the standard deviation for residuals (for normal errors only)
	# scale - the scale paramater for residuals (for cauchy errors only)
	# --------------------------------------------------------------------
	# OUTPUT: 
  	# noise - the simulated noisy responses
	# --------------------------------------------------------------------
	
	if (res.dist == "normal" & !is.na(sd)) {
		noisy <- f(x) + rnorm(n = length(x), mean = 0, sd = sd)
	} else if (res.dist == "cauchy" & !is.na(scale)) {
		noisy <- f(x) + rcauchy(n = length(x), location = 0, scale = scale)
	} else {
		stop(paste("Must give res.dist argument as ", 
			  "either \"normal\" or \"cauchy\", AND ",
			  "also specify sd (for normal) or ", 
			  "scale (for cauchy)", sep = ""))
	}
	
	return(noisy)
}

# ===========================================================================
# Cross-validation ==========================================================
# ===========================================================================

cv <- function(x.tr, y.tr, x.te, y.te, KERN.FUN, h) {
	# --------------------------------------------------------------------
	# Give cross validation mean squared prediction error 
	# --------------------------------------------------------------------
	# INPUTS:
	# x.tr - vector of predictors in * training * set
	# y.tr - vector of responses in * training * set
	# x.te - vector of predictors in * testing * set
	# y.te - vector of responses in * testing * set
	# KERN.FUN - some kernel function (e.g. Gaussian) 
	#            *** takes 2 arguments: distance (dist) and bandwidth (h)
	# h - vector or bandwidths
	# --------------------------------------------------------------------
	# OUTPUT: 
  	# mse - mean square predictive error
	# --------------------------------------------------------------------
	
	numbands <- length(h)
	mse <- rep(NA, numbands)
	
	for (i in 1:numbands) {
		y.pr <- sapply(
			x.te, 
			lin.smooth, 
			x = x.tr, 
			y = y.tr, 
			kern.fun = KERN.FUN, 
			h = h[i]
			)
			
		mse[i] <- mean((y.te - y.pr)^2)
	}
	
	return(mse)
}

# ===========================================================================
# Local polynomial regression ===============================================
# ===========================================================================

loc.lin <- function(x.new, x.vec, y.vec, h) {
	# --------------------------------------------------------------------
	# Give local linear estimator at some new point with normal kernel
	# --------------------------------------------------------------------
	# INPUTS:
	# x.new is some new point on x
	# x.vec - vector of x in sample
	# y.vec - vector of y in sample
	# h - bandwidth
	# --------------------------------------------------------------------
	# OUTPUT: 
  	# fit - the estimated value of f at x using local linear estimator
	# --------------------------------------------------------------------	
	
	kern.x <- kern.norm(x.new - x.vec, h)
	
	s1 <- sum(kern.x * (x - x.new))
	s2 <- sum(kern.x * (x - x.new)^2)
	# w.x <- kern.x * (s2 * x.new - s1 * (x - x.new))
	w.x <- kern.x * (s2 - s1 * (x - x.new))
	
	fit <- crossprod(w.x, y.vec) / sum(w.x)
	
	return(fit)
}

loc.pol <- function(x.new, x.vec, y.vec, D, h, give.mat = FALSE) {
	# --------------------------------------------------------------------
	# Local polynomial reg. estimator at some new point with normal kernel
	# --------------------------------------------------------------------
	# INPUTS:
	# x.new is some new point on x
	# x.vec - vector of x in sample
	# y.vec - vector of y in sample
	# D - dimension of the polynomial estimator (e.g., D = 1 means linear)
	# h - bandwidth for normal kernel
	# --------------------------------------------------------------------
	# OUTPUT (list): 
  	# fit - the estimated value of f at x using local linear estimator
	# hatmat.vec - normalized weights, to be used in creating a hat matrix
	# --------------------------------------------------------------------	
	
	# Number of observations
	N <- length(x.vec)
	
	# Vector of weights from Gaussian kernel
	w <- kern.norm(x.vec - x.new, h) / h
	w <- w / sum(w)
	
	# Create (non-normalized weights)
	if (D == 0) { # case of local constant estimator
		
		b <- matrix(w, nrow = 1)
		
	} else { # case of local polynomial estimator
		
		# Create R matrix
		R.x <- matrix(nrow = N, ncol = (D + 1))
		R.x[, 1] <- rep(1, N)
		
		for (j in 2:(D+1)) {
			R.x[, j] <- (x.vec - x.new)^(j - 1)
		} 
		
		W.diag <- diag(w)
		
		RxTW <- crossprod(R.x, W.diag)
		
		B <- solve(RxTW %*% R.x) %*% RxTW
		b <- matrix(B[1, ], nrow = 1)
	}
	
	fit <- tcrossprod( (b / sum(b)) , y.vec)
	
	# Should the hat matrix vector be given? 
	if (give.mat) {
		output <- list("fit" = fit, "hatmat.vec" = as.vector(b / sum(b)))
	} else {
		output <- fit
	}
	
	return(output)
}

loc.pol.hatmat <- function(x, y, D, h) {
	# --------------------------------------------------------------------
	# Create a hat matrix from local polynomial estimator
	# --------------------------------------------------------------------
	# INPUTS:
	# x - vector of x in sample
	# y - vector of y in sample
	# D - dimension of the polynomial estimator (e.g., D = 1 means linear)
	# h - bandwidth for normal kernel
	# --------------------------------------------------------------------
	# OUTPUT (list): 
  	# hatmat - hat matrix from local polynomial estimator
	# --------------------------------------------------------------------	
	
	N <- length(x)
	
	hatmat <- matrix(nrow = N, ncol = N)
	
	for (i in 1:N) {
		hatmat[i, ] <- loc.pol(x[i], x, y, D, h, give.mat = TRUE)$hatmat.vec
	}
	
	return(hatmat)
}

loocv <- function(y, H) {
	# --------------------------------------------------------------------
	# Generic leave-one-out cross validation
	# --------------------------------------------------------------------
	# INPUTS:
	# y - the response vector
	# H - the "hat" matrix
	# --------------------------------------------------------------------
	# OUTPUT: 
  	# loocv - leave-one-out cross validation error
	# --------------------------------------------------------------------	
	
	y.hat <- H %*% y
	
	loocv <- sum( ( {y - y.hat} / {1 - diag(H)} )^2 )
	
	return(loocv)
}

# ===========================================================================
# Gaussian process ==========================================================
# ===========================================================================

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
	
	# dimension of MVN
	p <- length(mu)
	
	# Check if inputs are valid (dimensions match, Sigma is square and p.s.d.)
	cond<- (ncol(Sigma) != p) | 
		   (nrow(Sigma) != p) | 
		   (max(eigen(Sigma)$values) <= 0) 
		   
	if (cond) {
		return("Try again...")
	}
	
	# Generate n*p univariate standard normal variables
	z     <- matrix(rnorm(n*p), nrow = p)
	
	# Create a matrix containing copies of mu
	mumat <- matrix(rep(mu, n), nrow = p)
	
	# Decompose Sigma into Sigma = L %*% Lt
	Lt <- chol(Sigma)
	
	# Generate sample with affine transformation of z
	x <- crossprod(Lt, z) + mumat
	
	return(t(x))
}

ell2 <- function(x) {
	# Compute the ell2 norm of x, a vector in Euclidean space
	
	return(sqrt(sum(x^2)))
}

C.SE <- function(dist, bool = F, params = NA) {
	# ----------------------------------------------------------------------
	#  Compute the (i, j) element of a squared exp. covariance matrix
	#  
	# ----------------------------------------------------------------------
	# INPUTS:
	# dist is some measure of distance between two points
	# params should be a vector of three hyperparameters
	#      1) b
	#      2) tau1.sq
	#      3) tau2.sq
	#
	# ----------------------------------------------------------------------
	# OUTPUT: 
  	# c.se is the value of the Matern-5/2 covariance matrix for x.i and x.j
	# ----------------------------------------------------------------------
	
	if (prod(is.na(params))) {
		return("Must have three valid parameters.")
	}
	
	if (length(params) != 3) {
		return("Must have three valid parameters.")
	}
	
	d <- dist
	
	b       <- params[1]
	tau1.sq <- params[2]
	tau2.sq <- params[3]
	
	# Euclidean distance between x.i and x.j
	# d <- ell2(x.i - x.j)
	
	c.se <- tau1.sq * exp(-0.5 * (d / b)^2) + tau2.sq * (bool)
	
	return(c.se)
}

# C.SE <- function(x.i, x.j, params = NA) {
# 	# ----------------------------------------------------------------------
# 	#  Compute the (i, j) element of a squared exp. covariance matrix
# 	#
# 	# ----------------------------------------------------------------------
# 	# INPUTS:
# 	# x.i and x.j are two vectors in same space (need not be [0, 1])
# 	# params should be a vector of three hyperparameters
# 	#      1) b
# 	#      2) tau1.sq
# 	#      3) tau2.sq
# 	#
# 	# ----------------------------------------------------------------------
# 	# OUTPUT:
#   	# c.se is the value of the Matern-5/2 covariance matrix for x.i and x.j
# 	# ----------------------------------------------------------------------
#
# 	if (prod(is.na(params))) {
# 		return("Must have three valid parameters.")
# 	}
#
# 	if (length(params) != 3) {
# 		return("Must have three valid parameters.")
# 	}
#
# 	b       <- params[1]
# 	tau1.sq <- params[2]
# 	tau2.sq <- params[3]
#
#
# 	b       <- params[1]
# 	tau1.sq <- params[2]
# 	tau2.sq <- params[3]
#
# 	# Euclidean distance between x.i and x.j
# 	d <- ell2(x.i - x.j)
#
# 	c.se <- tau1.sq * exp(-0.5 * (d / b)^2) + tau2.sq * (x.i == x.j)
#
# 	return(c.se)
# }


C.M52 <- function(dist, bool = F, params = NA) {
	# ----------------------------------------------------------------------
	#  Compute the (i, j) element of a Matern-5/2 covariance matrix
	#  
	# ----------------------------------------------------------------------
	# INPUTS:
	# dist is some measure of distance between two points
	# params should be a vector of three hyperparameters
	#      1) b
	#      2) tau1.sq
	#      3) tau2.sq
	#
	# ----------------------------------------------------------------------
	# OUTPUT: 
  	# c.m52 is the value of the Matern-5/2 covariance matrix for x.i and x.j
	# ----------------------------------------------------------------------
	
	if (prod(is.na(params))) {
		return("Must have three valid parameters.")
	}
	
	if (length(params) != 3) {
		return("Must have three valid parameters.")
	}
	
	d <- dist
	
	b       <- params[1]
	tau1.sq <- params[2]
	tau2.sq <- params[3]
	
	# Euclidean distance between x.i and x.j
	# d <- ell2(x.i - x.j)
	
	c.m52 <- tau1.sq * ( 1 + (5^0.5 * d / b) + (5 / 3 * (d / b)^2) ) * 
			 exp(-5^0.5 * d / b)  + tau2.sq * (bool)
	
	return(c.m52)
}

# C.M52 <- function(x.i, x.j, params = NA) {
# 	# ----------------------------------------------------------------------
# 	#  Compute the (i, j) element of a Matern-5/2 covariance matrix
# 	#
# 	# ----------------------------------------------------------------------
# 	# INPUTS:
# 	# x.i and x.j are two vectors in same space (need not be [0, 1])
# 	# params should be a vector of three hyperparameters
# 	#      1) b
# 	#      2) tau1.sq
# 	#      3) tau2.sq
# 	#
# 	# ----------------------------------------------------------------------
# 	# OUTPUT:
#   	# c.m52 is the value of the Matern-5/2 covariance matrix for x.i and x.j
# 	# ----------------------------------------------------------------------
#
# 	if (prod(is.na(params))) {
# 		return("Must have three valid parameters.")
# 	}
#
# 	if (length(params) != 3) {
# 		return("Must have three valid parameters.")
# 	}
#
# 	b       <- params[1]
# 	tau1.sq <- params[2]
# 	tau2.sq <- params[3]
#
# 	# Euclidean distance between x.i and x.j
# 	d <- ell2(x.i - x.j)
#
# 	c.m52 <- tau1.sq * ( 1 + (5^0.5 * d / b) + (5 / 3 * (d / b)^2) ) *
# 			 exp(-5^0.5 * d / b)  + tau2.sq * (x.i == x.j)
#
# 	return(c.m52)
# }

make.covmat <- function(x, cov.fun, params = NA) {
	# ----------------------------------------------------------------------
	#  Compute the symmetric covariance matrix for a GP, with some cov. function
	# ----------------------------------------------------------------------
	# INPUTS:
	# x is a vector of N values in [0, 1]
	# params should be a vector of three hyperparameters
	#      1) b
	#      2) tau1.sq
	#      3) tau2.sq
	#
	# ----------------------------------------------------------------------
	# OUTPUT: 
  	# covmat is the covariance matrix of GP
	# ----------------------------------------------------------------------
	
	if (prod(is.na(params))) {
		return("Must have three valid parameters.")
	}
	
	if (length(params) != 3) {
		return("Must have three valid parameters.")
	}
	
	N <- length(x)
	
	# covmat <- matrix(nrow = N, ncol = N)
	
	X <- matrix(rep(x, each = N), nrow = N)
	
	dist <- abs(X - t(X))
	
	covmat <- cov.fun(dist, bool = diag(N), params = params)
	
	# for (j in 1:N) {
	# 	for (i in j:N) {
	# 		dist <- ell2(x[i] - x[j])
	# 		covmat[i, j] <- cov.fun(dist,bool = (x[i] == x[j]), params = params)
	# 		covmat[j, i] <- covmat[i, j]
	# 	}
	# }
	
	return(covmat)
}


make.covmat2 <- function(x, cov.fun, params = NA, dist.params = NA) {
	# ----------------------------------------------------------------------
	#  Compute the symmetric covariance matrix for a GP, with some cov. function
	#  for x in R-2
	# ----------------------------------------------------------------------
	# INPUTS:
	# x - vector of x in sample 
	# cov.fun - covariance function
	# params - vector of three hyperparameters
	#      1) b
	#      2) tau1.sq
	#      3) tau2.sq
	# dist.params - a vector of two parameters for non-isotropic distance
	#
	# ----------------------------------------------------------------------
	# OUTPUT: 
  	# covmat is the covariance matrix of GP
	# ----------------------------------------------------------------------
	
	if (prod(is.na(params))) {
		return("Must have three valid parameters.")
	}
	
	if (length(params) != 3) {
		return("Must have three valid parameters.")
	}
	
	N <- nrow(x)
	
	X <- matrix(rep(x[, 1], each = N), nrow = N)
	Y <- matrix(rep(x[, 2], each = N), nrow = N)
	
	dist <- sqrt( (X - t(X))^2 / dist.params[1]^2 + 
				  (Y - t(Y))^2 / dist.params[2]^2 )
	
	
	covmat <- cov.fun(dist, bool = diag(N), params = params)
	
	return(covmat)
}

GP.predict <- function(x, y, x.new, cov.fun, params = NA, res.var) {
	# --------------------------------------------------------------------
	# Local polynomial reg. estimator at some new point with normal kernel
	# --------------------------------------------------------------------
	# INPUTS:
	# x - vector of x in sample
	# y - vector of y in sample
	# x.new - vector of y in sample
	# D - dimension of the polynomial estimator (e.g., D = 1 means linear)
	# h - bandwidth for normal kernel
	# --------------------------------------------------------------------
	# OUTPUT (list): 
  	# fit - the estimated value of f at x using local linear estimator
	# hatmat.vec - normalized weights, to be used in creating a hat matrix
	# --------------------------------------------------------------------	
	
	N <- length(x)
	M <- length(x.new)
	
	big.C <- make.covmat(c(x, x.new), cov.fun, params = params)
	
	C <- big.C[1:N, 1:N]
	C.s <- big.C[(N+1):(N+M), 1:N]
	C.ss <- big.C[(N+1):(N+M), (N+1):(N+M)]
	
	C.inv <- solve(C + res.var * diag(N))
	
	post.mean <- as.vector(C.s %*% C.inv %*% y)
	post.covm <- C.ss - C.s %*% C.inv %*% t(C.s)
	
	post.var <- diag(post.covm)
	
	mylist <- list("post.mean" = post.mean, "post.var"  = post.var)
	
	return(mylist)
}

GP.predict2 <- function(x, y, x.new, cov.fun, params = NA, 
	dist.params = c(1, 1), res.var) {
	# --------------------------------------------------------------------
	# Local polynomial reg. estimator at some new point with normal kernel
	# --------------------------------------------------------------------
	# INPUTS:
	# x - vector of x in sample
	# y - vector of y in sample
	# x.new - sequence of new points
	# params should be a vector of three hyperparameters
	#      1) b
	#      2) tau1.sq
	#      3) tau2.sq
	# dist.params - a vector of two parameters for non-isotropic distance
	# res.var - variance of residuals
	# --------------------------------------------------------------------
	# OUTPUT (list): 
  	# post.mean - posterior mean for f(x.new)
	# post.var - posterior covariance matrix for f(x.new)
	# --------------------------------------------------------------------	
	
	N <- nrow(x)
	M <- nrow(x.new)
	
	big.C <- make.covmat2(rbind(x, x.new), cov.fun, 
						  params = params, dist.params = dist.params)
	
	C <- big.C[1:N, 1:N]
	C.s <- big.C[(N+1):(N+M), 1:N]
	C.ss <- big.C[(N+1):(N+M), (N+1):(N+M)]
	
	C.inv <- solve(C + res.var * diag(N))
	
	post.mean <- as.vector(C.s %*% C.inv %*% y)
	post.covm <- C.ss - C.s %*% C.inv %*% t(C.s)
	
	post.var <- diag(post.covm)
	
	mylist <- list("post.mean" = post.mean, "post.var"  = post.var)
	
	return(mylist)
}
