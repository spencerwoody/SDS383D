library(mosaic)
library(MASS)
library(mvtnorm)
data(SaratogaHouses)  # ?SaratogaHouses for info

# for plotting
ghost_grey = rgb(0.5,0.5,0.5,0.2)

# Just pick some obvious variables
lm1 = lm(price ~ livingArea + lotSize + fuel + bedrooms + bathrooms, data=SaratogaHouses)
summary(lm1)

# note: mosaic overloads plot with an lm-like syntax that I prefer
plot(price ~ fitted(lm1), data=SaratogaHouses, pch=19, col=ghost_grey)
abline(0,1)

hist(resid(lm1))

# 4th (influence) plot: house 962 is a super high leverage point
plot(lm1)

SaratogaHouses[962,]
fitted(lm1)[962]
which.max(cooks.distance(lm1))

# Compare with rlm which uses Huber errors
lm2 = MASS::rlm(price ~ livingArea + lotSize + fuel + bedrooms + bathrooms, data=SaratogaHouses)
summary(lm2)

# some big differences in the fitted coefficients
cbind(coef(lm1), coef(lm2))

# and the fitted values
fitted(lm1)[962]
fitted(lm2)[962]


# Time for a Bayesian heavy-tailed error model

# Could extract the model matrix to go into the Bayesian linear model
# word to the wise: scale your covariates if you're
# going to use a common variance parameter for all betas...
X_raw = model.matrix(lm1)[,-1]
X = cbind(1, scale(X_raw))
head(X)

y <- SaratogaHouses$price

n <- nrow(X)
p <- ncol(X)

# -------------------------------------------------------------------------
# Specify prior model -----------------------------------------------------
# -------------------------------------------------------------------------

# Hyperparameters: K = diag(k1, k2), Lambda, d, eta, m

m <- matrix(rep(0, p), ncol = 1)
K <- diag(0.01, p)

Lambda <- diag(n)

d <- 0.02
eta <- 0.02

# For heavy tails only
h <- 2


# Fit non-heavy tail regression

XtLambda <- crossprod(X, Lambda)

K.star   <- crossprod(t(XtLambda), X) + K
m.star   <- solve(K.star, crossprod(t(XtLambda), y) + crossprod(K, m))
eta.star <- eta + crossprod(Lambda %*% y, y) + 
			crossprod(K %*% m, m) - crossprod(K.star %*% m.star, m.star)
d.star   <- d + n

m.star

# Prep for Gibbs sampler for model with heavy tails

nruns <- 5e3
burn  <- 1e3

beta.post <- matrix(nrow = nruns, ncol = p)
omega.post <- rep(NA, nruns)
lambda.post <- matrix(nrow = nruns, ncol = n)

beta.post[1, ] <- rep(0, p)
omega.post[1] <- 1
lambda.post[1, ] <- rep(1, n)

# Gibbs sampler
for (i in 2:nruns) {
	# Create lambda matrix
	Lambda.i <- diag(lambda.post[i - 1, ])
	
	# Precache this crossproduct
	XtLambda.i <- crossprod(X, Lambda.i)

	# Update parameters of posterior conditionals
	K.star.i   <- crossprod(t(XtLambda.i), X) + K
	m.star.i   <- solve(K.star.i, crossprod(t(XtLambda.i), y) + crossprod(K, m))
	eta.star.i <- eta + crossprod(Lambda.i %*% y, y) + 
			crossprod(K %*% m, m) - crossprod(K.star.i %*% m.star.i, m.star.i)
	
	# Draw from conditional posterior
	omega.post[i] <- rgamma(1, d.star / 2, eta.star.i / 2)
	beta.post[i, ] <- rmvnorm(1, m.star.i, 
					         sigma = solve(omega.post[i] * K.star.i))
	lambda.post[i, ] <- rgamma(n, rep(h/2 + 1/2, n), 
		h/2 + omega.post[i] * ( y- X %*% beta.post[i, ] )^2 / 2)
	
	print(i)
}

# Burn-out
beta.post <- beta.post[-(1:burn), ]
omega.post <- omega.post[-(1:burn)]
lambda.post <- lambda.post[-(1:burn), ] 

lambdas <- colMeans(lambda.post)

# Plot the lambdas
plot(sort(1 / lambdas), ylab = "1 / lambda")

# Mean posterior for betas
m.star2 <- matrix(colMeans(beta.post), ncol = 1)
rownames(m.star2) <- rownames(m.star)

# Frequentist

freq.lm <- my.lm(X, y)
freq.beta <- freq.lm$Beta.hat

freq.est <- matrix(freq.beta, ncol = 1)
freq.est

# Non-heavy tail Bayes
m.star 

# Heavy tail Bayes
m.star2

# Compare these estimates
comp <- cbind(freq.est, m.star, m.star2)
rownames(comp) <- rownames(m.star)
comp

plot(X %*% m.star, y)
abline(0, 1)

hist(y - X %*% m.star)

plot(X %*% m.star2, y)
abline(0, 1)

hist(y - X %*% m.star2)
