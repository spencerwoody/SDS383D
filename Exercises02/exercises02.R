###########################################################
######### Created by Spencer Woody on 04 Feb 2017 #########
###########################################################

library(ggplot2)
library(mvtnorm)
library(wesanderson)
library(extrafont)

source("myfuns02.R")

# Prep color pallette
pal <- wes_palette("Zissou", 5)
col1 <- pal[1]
col2 <- pal[4]
col3 <- pal[5]

# -------------------------------------------------------------------------
# Read in the data --------------------------------------------------------
# -------------------------------------------------------------------------

GDP <- read.csv("gdpgrowth.csv", header = T)

y <- matrix(GDP$GR6096, ncol = 1)
X <- matrix(GDP$DEF60, ncol = 1)

n <- length(y)

X <- cbind(rep(1, n), X)
colnames(X) <- c("int", "DEF60")

# -------------------------------------------------------------------------
# Specify prior model -----------------------------------------------------
# -------------------------------------------------------------------------

freq.lm <- my.lm(X, y)
freq.beta <- freq.lm$Beta.hat

# -------------------------------------------------------------------------
# Specify prior model -----------------------------------------------------
# -------------------------------------------------------------------------

# Hyperparameters: K = diag(k1, k2), Lambda, d, eta, m

k1 <- 0.01
k2 <- 0.01

m <- matrix(c(0, 0), ncol = 1)
K <- diag(c(k1, k2))

Lambda <- diag(n)

d <- 0.02
eta <- 0.02

# For heavy tails only
h <- 0.02

# -------------------------------------------------------------------------
# Obtain posterior model (no heavy tails) ---------------------------------
# -------------------------------------------------------------------------

XtLambda <- crossprod(X, Lambda)

K.star   <- crossprod(t(XtLambda), X) + K
m.star   <- solve(K.star, crossprod(t(XtLambda), y) + crossprod(K, m))
eta.star <- eta + crossprod(Lambda %*% y, y) + 
			crossprod(K %*% m, m) - crossprod(K.star %*% m.star, m.star)
d.star   <- d + n

m.star

# Prep for Gibbs sampler for model with heavy tails

nruns <- 1e4
burn  <- 2e3

beta.post <- matrix(nrow = nruns, ncol = 2)
omega.post <- rep(NA, nruns)
lambda.post <- matrix(nrow = nruns, ncol = n)

beta.post[1, ] <- c(0, 0)
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
}

# Burn-out
beta.post <- beta.post[-(1:burn), ]
omega.post <- omega.post[-(1:burn)]
lambda.post <- lambda.post[-(1:burn), ] 

# Mean posterior for betas
m.star2 <- colMeans(beta.post)

# Plot Frequentist, Bayesian, and Bayesian heavy-tail linear models
r <- ggplot(GDP, aes(DEF60, GR6096)) + geom_point() +
     geom_abline(mapping = aes(colour = "Frequentist", 
	 intercept = freq.beta[1], 
	 slope = freq.beta[2]), 
	 show.legend = T) + 
	 geom_abline(mapping = aes(colour = "Bayesian heavy tails", 
	 intercept = m.star2[1], 
	 slope = m.star2[2]), 
	 show.legend = T) +
	 geom_abline(mapping = aes(colour = "Bayesian iid residuals", 
	 intercept = m.star[1], 
	 slope = m.star[2]), 
	 show.legend = T) +
	 xlab("Defense spending as fraction of GDP (DEF60)") +
	 ylab("GDP Growth Rate (GR6096)") +
	 ggtitle("Defense spending vs. GDP Growth Rate") + 
	 scale_color_manual(name = "Method",
	 values = c("Frequentist" = col1, 
	 "Bayesian iid residuals" = col2, "Bayesian heavy tails" = col3)) +
	 theme(legend.position = c(0.85, 0.15),
	 text = element_text(family="Palatino"))
r
	 
# Save this plotr
pdf("img/compplot.pdf", width = 8, height = 6)	
r
dev.off()

# Plot of reciprocals of lambdas
plot(sort(1 / colMeans(lambda.post)))

plot(sort(colMeans(lambda.post)))

# Plot gibbs sampler results
plot(omega.post,
	type = "l")

q <- qplot(1:length(beta.post[, 1]), beta.post[, 1], geom = "blank") + 
	 geom_line() + 
	 ggtitle("Traceplot of Gibbs sampler for intercept term") +
	 xlab("iteration") +
	 ylab("Intercept") + 
	 theme(text = element_text(family="Palatino"))

w <- qplot(1:length(beta.post[, 2]), beta.post[, 1], geom = "blank") + 
	 geom_line() + 
	 ggtitle("Traceplot of Gibbs sampler for DEF60 term") +
	 xlab("iteration") +
	 ylab("Intercept") + 
	 theme(text = element_text(family="Palatino"))


