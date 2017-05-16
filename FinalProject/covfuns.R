
C.M52 <- function(dist, bool = F, params = NA) {
	# ----------------------------------------------------------------------
	#  Compute the (i, j) element of a Matern(5/2) covariance matrix
	#  k(t, t') = tau1.sq * exp(1) * exp tau2.sq * (t == t')
	# ----------------------------------------------------------------------
	# INPUTS:
	# dist - Euclidean distance between two points
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
	
	c.m52 <- tau1.sq * { 1 + (5^0.5 * d / b) + (5 / 3 * (d / b)^2) } * 
			 exp(-5^0.5 * d / b)  + tau2.sq * (bool)
	
	return(c.m52)
}

make.covmat <- function(t1, t2 = NA, cov.fun, params = NA) {
	# ----------------------------------------------------------------------
	#  Compute the Covariance matrix for a GP, with some cov. function
	# ----------------------------------------------------------------------
	# INPUTS:
	# t1 - a vector of times
	# t2 - a vector of times
	# cov.fun - some covariance function
	# params - a vector of three hyperparameters
	#      1) b
	#      2) tau1.sq
	#      3) tau2.sq
	#
	# ----------------------------------------------------------------------
	# OUTPUT: 
  	# covmat is the covariance matrix of GP
	# ----------------------------------------------------------------------
	
	if (sum(is.na(t2))) {
		t2 <- t1
	} 
	
	
	# N <- length(x)
	
	T1 <- matrix(rep(t1, length(t2)), byrow = F, nrow = length(t1))
	T2 <- matrix(rep(t2, each = length(t1), byrow = T), ncol = length(t2))
	
	# covmat <- matrix(nrow = N, ncol = N)
	
	dist <- abs(T1 - T2)
	
	bool <- (T1 == T2) * 1
	
	covmat <- cov.fun(dist, bool = bool, params = params)
	
	return(covmat)
}

neg.ll <- function(params, psi.n, z.n) {
	# ----------------------------------------------------------------------
	#  Compute the negative log-likelihood of psi-vector 
	# ----------------------------------------------------------------------
	theta.g   <- c(exp(params[1]), exp(params[2]), 0)
	theta.p   <- c(exp(params[3]), exp(params[3]), 0)
	
	K.g <- make.covmat(t.n, cov.fun = C.M52, params = theta.g)
	K.p <- make.covmat(t.nr, cov.fun = C.M52, params = theta.p)
	K.n <- K.g + bdiag(rep(list(K.p), 3))

	K.nINV <- solve(K.n)

	Sigma.n <- as.matrix(solve(K.nINV + Omega.n))
	
	out <- - dmvnorm(
		psi.n, 
		mean = Sigma.n %*% z.n, 
		sigma = Sigma.n,
		log = TRUE)
	
	return(out)
}

