###########################################################
######### Created by Spencer Woody on 27 Mar 2017 #########
###########################################################

library(ggplot2)
library(RColorBrewer) # display.brewer.all()
library(Matrix)
library(sparseMVN)
library(mvtnorm)
library(gridExtra)

### --------------------------------------------------------------------------
### Data prep
### --------------------------------------------------------------------------

polls <- read.csv("Polls/polls.csv", header = T)

states <- polls$state

N <- nrow(polls)

X <- cbind(rep(1, n),
polls$edu == "HS", polls$edu == "SomeColl", polls$edu == "Bacc",
polls$age == "30to44", polls$age == "45to64", polls$age == "65plus",
polls$female,
polls$black)

y <- polls$bush

colnames(X) <- c("(Intercept)", "HS", "SomeColl", "Bacc", "30to44", "45to64",
"65plus", "female", "black")

W <- X

I <- length(unique(states))

p <- ncol(X)
q <- ncol(W)

X.list <- list()
W.list <- list()

for (i in unique(states)) {
	indices <- which(states == i)
	X.list[[ length(X.list) + 1 ]] <- X[indices, ]
	W.list[[ length(W.list) + 1 ]] <- W[indices, ]
}

X.block <- bdiag(X.list)
Xt.block <- t(X.block)
XtX.block <- crossprod(X.block)

W.block <- bdiag(X.list)
Wt.block <- t(W.block)
WtW.block <- crossprod(W.block)

SSX <- matrix(0, nrow = p, ncol = p) 

for (i in 1:I) {
	SSX <- SSX + XtX.block[((i - 1) * p + 1):(i * p), ((i - 1) * p + 1):(i * p)]
	# SSX <- SSX + crossprod(X[inds == i, ])
}

SSX.inv <- as.matrix(solve(SSX))





which.0 <- which(y == 0)
which.1 <- which(y == 1)
n0 <- length(which.0)
n1 <- length(which.1)

### --------------------------------------------------------------------------
### Gibbs sampler
### --------------------------------------------------------------------------

# Number of iterations and burn-in
niter <- 6000
nburn <- 1000

# Prep prepare to store samples
b <- matrix(nrow = niter, ncol = I * q)
beta <- matrix(nrow = niter, ncol = p)
D <- array(dim = c(p, p, niter))
# z <- matrix(nrow = niter, ncol = N)

y.predmat <- matrix(nrow = niter - 1, ncol = N)

# Prior parameters for D
nu <- q + 1
Psi <- diag(q)

# Initialize the Gibbs sampler
b[1, ] <- rep(0.1, q * I)
beta[1, ] <- rep(0.1, p)
D[, , 1] <- diag(q)

# Initialize latent variables
z.n <- rnorm(N)

# Run the Gibbs sampler
for (iter in 2:niter) {
	# Change "current" values
	b.c <- b[iter - 1, ]
	beta.c <- beta[iter - 1, ]
	D.c <- D[, , iter - 1]
	z.c <- z.n
	
	# Update b
	u <- z.c - X %*% beta.c
	
	b.prec <- WtW.block + bdiag(rep(list(D.c), I))
	b.mean <- solve(b.prec, crossprod(W.block, u))

	mychol <- Cholesky(b.prec)
	
	b.n <- as.numeric(rmvn.sparse(1, b.mean, mychol, TRUE))
	
	# Update D
	B.n <- matrix(as.numeric(b.n), nrow = q)
	D.n <- rWishart(1, nu + I, Psi + tcrossprod(B.n))
	
	# Update beta
	v <- z.c - W.block %*% b.n
	beta.n <- as.numeric(rmvnorm(1, SSX.inv %*% crossprod(X, v), SSX.inv))
	
	# Generate z
	z.meanvec <- as.numeric(X %*% beta.n + W.block %*% b.n)
	z.pred <- rnorm(N, z.meanvec, rep(1, N))
	y.predmat[iter - 1, ] <- (z.pred > 0) * 1
	
	z.n <- rep(0, N)
	z.n[which.0] <- rtnorm(n0, b = 0, mean = z.meanvec[which.0],sd = rep(1, n0))
	z.n[which.1] <- rtnorm(n1, a = 0, mean = z.meanvec[which.1],sd = rep(1, n1))
	
	# Progress report
	if (iter %% 100 == 0) {
		print(sprintf("Iteration %i out of %i...", iter, niter))
	}
	
	# Add updated parameter
	b[iter, ] <- b.n
	beta[iter, ] <- beta.n
	D[, , iter] <- D.n
}

# Remove burn-in
b <- b[-(1:nburn), ]
beta <- beta[-(1:nburn), ]
D <- D[, , -(1:nburn)]

length.post <- length(lambda)


plot(beta[, 1], type = "l")

# Confusion matrix
y.pred == 1 * y == 0


rtnorm <- function(n, a = -Inf, b = Inf, mean = 0, sd = 1) {
	
	Fa <- pnorm(a, mean, sd)
	Fb <- pnorm(b, mean, sd)
	
	u <- runif(n, 0, 1)
	x <- qnorm(u * (Fb - Fa) + Fa, mean, sd)
	
	return(x)
}





