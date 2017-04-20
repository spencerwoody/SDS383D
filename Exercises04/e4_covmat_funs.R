
C.SE <- function(dist, bool = F, params = NA) {
	# ----------------------------------------------------------------------
	#  Compute the (i, j) element of a squared exp. covariance matrix
	#  k(t, t') = tau1.sq * exp(- dist / b) + tau2.sq * (t == t')
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
	
	c.se <- tau1.sq * exp(- (d / b)^2 / 2) + tau2.sq * (bool)
	
	return(c.se)
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

full.covmat <- function(Yi, params) {
	# ----------------------------------------------------------------------
	#  Compute the full covariance matrix for gene expression data
	# ----------------------------------------------------------------------
	# INPUTS:
	# Yi - Dataframe from all replicates of all genes in the cluster, 
	#      (columns)
	#      1) group
	#      2) gene
	#      3) replicate
	# 	   4) time
	#      5) log2exp
	# params - a vector of seven parameters
	#      1) gamma.f - 
	#      2) alpha.f - 
	#      3) gamma.g - 
	#      4) alpha.g - 
	#      5) gamma.h - 
	#      6) alpha.h - 
	#      7) sigma2  - 
	#
	# ----------------------------------------------------------------------
	# OUTPUT: 
  	# negative log-likelihood 
	# ----------------------------------------------------------------------
	
	require(mvtnorm)
	
	# ----------------------------------------------------------------------
	# Renames the parameters
	# ----------------------------------------------------------------------
	
	gamma.f <- exp(params[1])
	alpha.f <- exp(params[2])
	gamma.g <- exp(params[3])
	alpha.g <- exp(params[4])
	gamma.h <- exp(params[5])
	alpha.h <- exp(params[6])
	sigma2  <- exp(params[7])
	
	
	# Number of genes in the cluster
	Ni <- length(unique(Yi$gene))

	# Choose one gene to make Sigma.n, which is same for all genes in cluster
	yn <- Yi[Yi$gene == Yi$gene[1], ]
	
	reps <- unique(Yi$replicate)
	numreps <- length(reps)
	
	# ----------------------------------------------------------------------
	# Make Kf block matrices
	# ----------------------------------------------------------------------
	
	tnr.list <- list()
	Kf.list <- list()

	for (k in 1:numreps) {
		tnr.list[[ k ]] <- yn$time[yn$replicate == reps[k]]
		Kf.list[[ k ]] <- make.covmat(
			tnr.list[[ k ]],
			cov.fun = C.SE,
			params = c(gamma.f, alpha.f, 0))
	}
	
	tn <- unlist(tnr.list) # Vector of times for all replicates in one gene
	D <- length(tn) # Number of observations across replicates for one gene
	
	
	# ----------------------------------------------------------------------
	# Make other parts of full likelihood covariance matrix
	# ----------------------------------------------------------------------
	
	# Within-cluster covariance
	Kh <- make.covmat(
		tn,
		cov.fun = C.SE,
		params = c(gamma.h, alpha.h, 0))

	Kh.long <- matrix(rep(Kh, Ni), nrow = nrow(Kh)) # make D x D blocks of Kh
	Kh.big <- do.call("rbind", rep(list(Kh.long), Ni))

	# Within-gene covariance
	Kg <- make.covmat(
		tn,
		cov.fun = C.SE,
		params = c(gamma.g, alpha.g, 0))

	# Full covariance matrix for one gene
	Sigma.n <- Kg + bdiag(Kf.list) + sigma2 * Diagonal(D)
	
	# Make covariance matrix for ALL genes in the cluster
	COVMAT <- as.matrix(Kh.big + bdiag(rep(list(Sigma.n), Ni)))
	
	return(COVMAT)
}

neg.ll <- function(Yi, params) {
	# ----------------------------------------------------------------------
	#  Compute the full covariance 
	# ----------------------------------------------------------------------
	COVMAT <- full.covmat(Yi, params)
	
	neg.ll <- - dmvnorm(Yi$log2exp, mean = rep(0, nrow(Yi)),
				 	  sigma = COVMAT, log = TRUE) 
					  
	return(neg.ll)
}

