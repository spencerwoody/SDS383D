###########################################################
######### Created by Spencer Woody on 27 Mar 2017 #########
###########################################################

rm(list=ls())

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
polls <- na.omit(polls)

states <- polls$state
states.list <- unique(states)

N <- nrow(polls)

# Design matrix with ordered data
X <- cbind(rep(1, N),
polls$edu %in% c("HS", "SomeColl", "Bacc"), 
polls$edu %in% c("SomeColl", "Bacc"), 
polls$edu %in% c("Bacc"),
polls$age %in% c("30to44", "45to64", "65plus"), 
polls$age %in% c("45to64", "65plus"), 
polls$age %in% c("65plus"),
polls$female,
polls$black)

# Response vector
y <- polls$bush

colnames(X) <- c("(Intercept)", "HS", "SomeColl", "Bacc", "30to44", "45to64",
"65plus", "female", "black")

W <- as.matrix(X[, 1])

I <- length(unique(states))

p <- ncol(X)
q <- ncol(W)

X.list <- list()
W.list <- list()

for (i in states.list) {
	indices <- which(states == i)
	X.list[[ length(X.list) + 1 ]] <- X[indices, ]
	W.list[[ length(W.list) + 1 ]] <- as.matrix(W[indices, ])
}

X.block <- bdiag(X.list)
Xt.block <- t(X.block)
XtX.block <- crossprod(X.block)

W.block <- bdiag(W.list)
Wt.block <- t(W.block)
WtW.block <- crossprod(W.block)

SSX <- matrix(0, nrow = p, ncol = p) 

for (i in 1:I) {
	SSX <- SSX + XtX.block[((i - 1) * p + 1):(i * p), ((i - 1) * p + 1):(i * p)]
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
niter <- 1000
nburn <- 100
length.post <- niter - nburn

# Prep prepare to store samples
b <- matrix(nrow = niter, ncol = I * q)
beta <- matrix(nrow = niter, ncol = p)
D <- array(dim = c(q, q, niter))

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
	# z.meanvec <- as.numeric(X %*% beta.n + W.block %*% b.n)
	# z.pred <- rnorm(N, z.meanvec, rep(1, N))
	# y.predmat[iter - 1, ] <- (z.pred > 0) * 1 # 
	z.pred <- as.numeric(X %*% beta.n + W.block %*% b.n)
	y.predmat[iter - 1, ] <- pnorm(z.pred) # possible fix
	
	z.n <- rep(NA, N)
	z.n[which.0] <- rtnorm(n0, b = 0, mean = z.pred[which.0],sd = rep(1, n0))
	z.n[which.1] <- rtnorm(n1, a = 0, mean = z.pred[which.1],sd = rep(1, n1))
	
	# Progress report
	if (iter %% 100 == 0) {
		print(sprintf("Iteration %i out of %i...", iter, niter))
	}
	
	# Add updated parameter
	b[iter, ] <- b.n
	beta[iter, ] <- beta.n
	D[, , iter] <- D.n
}

plot(beta[, 1], type = "l")

# Accuracy
y.pred <- colMeans(y.predmat)
sum((y == 1) * (y.pred > 0.5)) / N

# 95% Credible intervals for beta
betasummary <- t(as.matrix(apply(beta, 2, quantile, c(0.025, 0.500, 0.975))))
rownames(betasummary) <- colnames(X)
betasummary

# 95% Credible intervals for b
b.summary <- t(as.matrix(apply(b, 2, quantile, c(0.025, 0.500, 0.975))))

# Intercept
int.df <- data.frame(state = factor(states.list),
lo = b.summary[, 1],
md = b.summary[, 2],
hi = b.summary[, 3])

int.df <- int.df[order(int.df$md), ]
int.df$state <- factor(int.df$state, 
	levels = int.df$state)

intCI1 <- ggplot(int.df, aes(ymin = lo, ymax = hi, x = state)) + 
geom_linerange() + 
xlab("State") +
ylab("State-level Intercept") +
ggtitle("95% Credible intervals for state-level intercept") +
geom_hline(yintercept = 0, linetype = "dashed") +
geom_point(aes(y = md), col = "firebrick3") +
theme_bw() + 
theme(plot.title = element_text(hjust = 0.5),
text=element_text(family="Courier New", face="bold"))

intCI2 <- ggplot(int.df, aes(ymin = lo, ymax = hi, x = state)) + 
geom_linerange() + 
xlab("State") +
ylab("State-level Intercept") +
ggtitle("95% Credible intervals for state-level intercept") +
geom_hline(yintercept = 0, linetype = "dashed") +
geom_point(aes(y = md), col = "firebrick3") +
theme_bw() + 
theme(plot.title = element_text(hjust = 0.5))



rtnorm <- function(n, a = -Inf, b = Inf, mean = 0, sd = 1) {
	
	Fa <- pnorm(a, mean, sd)
	Fb <- pnorm(b, mean, sd)
	
	u <- runif(n, 0, 1)
	x <- qnorm(u * (Fb - Fa) + Fa, mean, sd)
	
	return(x)
}





