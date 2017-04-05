###########################################################
######### Created by Spencer Woody on 27 Mar 2017 #########
###########################################################

library(ggplot2)
library(RColorBrewer) # display.brewer.all()
library(sparseMVN)
library(mvtnorm)
library(gridExtra)

col0 <- "dodgerblue4"
col1 <- "darkorange"

### --------------------------------------------------------------------------
### Data prep
### --------------------------------------------------------------------------

# Read in the data
cheese.csv <- read.csv("Cheese/cheese.csv")

# Transform the data
cheese <- data.frame(store = factor(cheese.csv$store),
	logprice = log(cheese.csv$price),
	disp = cheese.csv$disp,
	logprice.disp = log(cheese.csv$price) * cheese.csv$disp,
	logvol = log(cheese.csv$vol))

# Save names of stores
stores.list <- unique(cheese$store)

# Change the levels of store to integers
levels(cheese$store) <- 1:length(levels(cheese$store))
cheese$store <- as.numeric(cheese$store)

### --------------------------------------------------------------------------
### EDA
### --------------------------------------------------------------------------

# OLS estimate & plot
mylm <- lm(logvol ~ logprice + disp + logprice.disp, data = cheese)
betahat <- coef(mylm)

OLSplot <- ggplot(cheese, aes(logprice, logvol)) + 
geom_abline(slope = betahat[2] + betahat[4], intercept = betahat[1] + betahat[3], col = col1) +
geom_abline(slope = betahat[2], intercept = betahat[1], col = "blue") + 
xlab("log-Price (USD)") +
ylab("log-Volume sold (units)") +
ggtitle("Demand curve for all stores with OLS line") +
geom_point(pch = 1, aes(colour = factor(disp))) +
scale_colour_manual("Display", values = c("blue", col1)) + 
theme(plot.title = element_text(hjust = 0.5)) 

pdf("Cheese/img/OLSplot.pdf")
OLSplot
dev.off()

OLSdf <- data.frame(
	store = factor(1:length(stores.list)), 
	intercept = as.numeric(rep(betahat[1], length(stores.list))),
	logprice = as.numeric(rep(betahat[2], length(stores.list))),
	disp = as.numeric(rep(betahat[3], length(stores.list))),
	logprice.disp = as.numeric(rep(betahat[4], length(stores.list)))
	)

facetplot <- ggplot(cheese, aes(logprice, logvol)) + 
geom_point(pch = 1, aes(colour = factor(disp))) + 
xlab("log-Price (USD)") +
ylab("log-Volume sold (units)") +
facet_wrap(~factor(store), ncol = 9) +
scale_colour_manual("Display", values = c("blue", col1))
	
OLSfacetplot <- facetplot + 
ggtitle("Demand curve for all stores with OLS line") +
geom_abline(data = OLSdf, aes(slope = logprice, intercept = intercept), col = "blue") +
geom_abline(data = OLSdf, aes(slope = logprice + logprice.disp, intercept = intercept + disp), col = col1) +
theme(# strip.background = element_blank(),
#        strip.text.x = element_blank(),
	   plot.title = element_text(hjust = 0.5))

pdf("Cheese/img/OLSfacetplot.pdf", width = 10, height = 10)
OLSfacetplot
dev.off()

### --------------------------------------------------------------------------
### Create matrices and block diagonal matrices
### --------------------------------------------------------------------------

cheese.mat <- as.matrix(cheese)

N <- nrow(cheese.mat) # Total number of observations

inds <- cheese.mat[, 1]
X <- cbind(rep(1, N), cheese.mat[, 2:4])
Z <- cbind(rep(1, N), cheese.mat[, 2:4])
y <- cheese.mat[, 5]

I <- length(unique(inds)) # Number of groups

p <- ncol(X)
q <- ncol(Z)

# Make block diagonal matrices
X.list <- list()
Z.list <- list()

for (i in unique(inds)) {
	X.list[[ length(X.list) + 1 ]] <- X[inds == i, ]
	Z.list[[ length(Z.list) + 1 ]] <- Z[inds == i, ]
}

X.block <- bdiag(X.list)
Xt.block <- t(X.block)
XtX.block  <- crossprod(X.block)

Z.block <- bdiag(Z.list)
Zt.block <- t(Z.block)
ZtZ.block <- crossprod(Z.block)

# Precache SSXinv
SSX <- matrix(0, nrow = p, ncol = p) 

for (i in 1:I) {
	SSX <- SSX + XtX.block[((i - 1) * p + 1):(i * p), ((i - 1) * p + 1):(i * p)]
	# SSX <- SSX + crossprod(X[inds == i, ])
}

SSX.inv <- as.matrix(solve(SSX))

### --------------------------------------------------------------------------
### Gibbs sampler
### --------------------------------------------------------------------------

# Number of iterations and burn-in
niter <- 6000
nburn <- 1000

# Prep prepare to store samples
b <- matrix(nrow = niter, ncol = I * q)
beta <- matrix(nrow = niter, ncol = p)
lambda <- rep(NA, niter)
D <- array(dim = c(p, p, niter))

# Prior parameters for D
nu <- q + 1
Psi <- diag(q)

# Initialize the Gibbs sampler
b[1, ] <- rep(0, q * I)
beta[1, ] <- rep(0, p)
D[, , 1] <- diag(q)
lambda[1] <- 1

# Run the Gibbs sampler
for (iter in 2:niter) {
	# Change "current" values
	b.c <- b[iter - 1, ]
	lambda.c <- lambda[iter - 1]
	beta.c <- beta[iter - 1, ]
	D.c <- D[, , iter - 1]
	
	# Update lambda
	RSS <- sum((y - X %*% beta.c - Z.block %*% b.c)^2)
	lambda.n <- rgamma(1, N / 2, RSS / 2)
	
	# Update b
	v <- y - X %*% beta.c
	
	b.prec <- lambda.n * ZtZ.block + bdiag(rep(list(D.c), I))
	b.mean <- lambda.n * solve(b.prec, crossprod(Z.block, v))

	mychol <- Cholesky(b.prec)
	
	b.n <- as.numeric(rmvn.sparse(1, b.mean, mychol, TRUE))
	
	# Update D
	B.n <- matrix(as.numeric(b.n), nrow = q)
	D.n <- rWishart(1, nu + I, Psi + tcrossprod(B.n))
	
	# Update beta
	w <- y - Z.block %*% b.n
	beta.n <- rmvnorm(1, SSX.inv %*% crossprod(X, w), SSX.inv / lambda.n)
	
	# Progress report
	if (iter %% 100 == 0) {
		print(sprintf("Iteration %i out of %i...", iter, niter))
	}
	
	# Add updated parameter
	b[iter, ] <- b.n
	lambda[iter] <- lambda.n
	beta[iter, ] <- beta.n
	D[, , iter] <- D.n
}

# Remove burn-in
b <- b[-(1:nburn), ]
lambda <- lambda[-(1:nburn)]
beta <- beta[-(1:nburn), ]
D <- D[, , -(1:nburn)]

length.post <- length(lambda)

# library(microbenchmark)
#
# microbenchmark(mychol <- Cholesky(b.prec))
# microbenchmark(b.mean <- lambda.n * solve(b.prec, crossprod(Z.block, v)))
# microbenchmark(covmat <- solve(b.prec))
# microbenchmark(samp <- rmvnorm(1, b.mean, as.matrix(covmat)))
# microbenchmark(D.n <- rWishart(1, nu + I, Psi + tcrossprod(B.n)))



### --------------------------------------------------------------------------
### Traceplots and histograms
### --------------------------------------------------------------------------

# Check mixing
lambdamix <- qplot(1:length.post,lambda, geom = "path") + 
xlab("iteration") +
ylab("lambda") +
ggtitle("Traceplot for Gibbs sampler for lambda") +
theme(plot.title = element_text(hjust = 0.5))
lambdamix
	
intmix <- qplot(1:length.post, beta[, 1], geom = "path") + 
xlab("iteration") +
ylab("intercept") +
ggtitle("Traceplot for Gibbs sampler for intercept") +
theme(plot.title = element_text(hjust = 0.5))
intmix

logpricemix <- qplot(1:length.post, beta[, 2], geom = "path") + 
xlab("iteration") +
ylab("logprice") +
ggtitle("Traceplot for Gibbs sampler for logprice") +
theme(plot.title = element_text(hjust = 0.5))
logpricemix

dispmix <- qplot(1:length.post, beta[, 3], geom = "path") + 
xlab("iteration") +
ylab("disp") +
ggtitle("Traceplot for Gibbs sampler for disp") +
theme(plot.title = element_text(hjust = 0.5))
dispmix

logprice.dispmix <- qplot(1:length.post, beta[, 4], geom = "path") + 
xlab("iteration") +
ylab("logprice.disp") +
ggtitle("Traceplot for Gibbs sampler for logprice.disp") +
theme(plot.title = element_text(hjust = 0.5))
logprice.dispmix

# Histograms

lambdahist <- qplot(lambda, geom = 'blank') +   
geom_histogram(aes(y = ..density..), 
    alpha = 0.8,  
	bins = 20, 
	colour = "grey95", 
	fill = "grey20") +
xlab("lambda") +
ylab("Density") +
ggtitle("Posterior draws for lambda") +
theme_bw() +
theme(plot.title = element_text(hjust = 0.5))
lambdahist

inthist <- qplot(beta[, 1], geom = 'blank') +   
geom_histogram(aes(y = ..density..), 
    alpha = 0.8,  
	bins = 20, 
	colour = "grey95", 
	fill = "grey20") +
xlab("intercept") +
ylab("Density") +
ggtitle("Posterior draws for intercept") +
theme_bw() +
theme(plot.title = element_text(hjust = 0.5))
inthist

logpricehist <- qplot(beta[, 2], geom = 'blank') +   
geom_histogram(aes(y = ..density..), 
    alpha = 0.8,  
	bins = 20, 
	colour = "grey95", 
	fill = "grey20") +
xlab("logprice") +
ylab("Density") +
ggtitle("Posterior draws for logprice") +
theme_bw() +
theme(plot.title = element_text(hjust = 0.5))
logpricehist

disphist <- qplot(beta[, 3], geom = 'blank') +   
geom_histogram(aes(y = ..density..), 
    alpha = 0.8,  
	bins = 20, 
	colour = "grey95", 
	fill = "grey20") +
xlab("disp") +
ylab("Density") +
ggtitle("Posterior draws for disp") +
theme_bw() +
theme(plot.title = element_text(hjust = 0.5))
disphist

logprice.disphist <- qplot(beta[, 4], geom = 'blank') +   
geom_histogram(aes(y = ..density..), 
    alpha = 0.8,  
	bins = 20, 
	colour = "grey95", 
	fill = "grey20") +
xlab("lambda") +
ylab("logprice.disp") +
ggtitle("Posterior draws for logprice.disp") +
theme_bw() +
theme(plot.title = element_text(hjust = 0.5))
logprice.disphist

### --------------------------------------------------------------------------
### Credible intervals ("grand mean")
### --------------------------------------------------------------------------

# 95% credible intervals
# beta
betasummary <- t(as.matrix(apply(beta, 2, quantile, c(0.025, 0.500, 0.975))))
rownames(betasummary) <- names(betahat)
betasummary

# lambda
lambdasummary <- quantile(lambda, c(0.025, 0.500, 0.975))
lambdasummary

### --------------------------------------------------------------------------
### Credible intervals plots for subject level effects
### --------------------------------------------------------------------------

# Subject level effects
bsummary <- t(as.matrix(apply(b, 2, quantile, c(0.025, 0.500, 0.975))))

int.summary <- bsummary[seq(1, 352, by = 4), ]
logprice.summary <- bsummary[seq(2, 352, by = 4), ]
disp.summary <- bsummary[seq(3, 352, by = 4), ]
logprice.disp.summary <- bsummary[seq(4, 352, by = 4), ]

# Intercept
int.summary.df <- data.frame(store = factor(1:88), 
lo = int.summary[, 1],
md = int.summary[, 2],
hi = int.summary[, 3])

int.summary.df <- int.summary.df[order(int.summary.df$md), ]
int.summary.df$store <- factor(int.summary.df$store, 
	levels = int.summary.df$store)

intCI <- ggplot(int.summary.df, aes(ymin = lo, ymax = hi, x = store)) + 
geom_linerange() + 
xlab("Store number") +
ylab("Store-level Intercept") +
ggtitle("95% Credible intervals for store-level intercept") +
geom_hline(yintercept = 0, linetype = "dashed") +
geom_point(aes(y = md), col = "firebrick3") +
theme_bw() + 
theme(plot.title = element_text(hjust = 0.5))

# logprice
logprice.summary.df <- data.frame(store = factor(1:88), 
lo = logprice.summary[, 1],
md = logprice.summary[, 2],
hi = logprice.summary[, 3])

logprice.summary.df <- logprice.summary.df[order(logprice.summary.df$md), ]
logprice.summary.df$store <- factor(logprice.summary.df$store, 
	levels = logprice.summary.df$store)

logpriceCI <- ggplot(logprice.summary.df,aes(ymin = lo, ymax = hi,x = store)) + 
geom_linerange() + 
xlab("Store number") +
ylab("Store-level logprice term") +
ggtitle("95% Credible intervals for store-level logprice term") +
geom_hline(yintercept = 0, linetype = "dashed") +
geom_point(aes(y = md), col = "firebrick3") +
theme_bw() + 
theme(plot.title = element_text(hjust = 0.5))
	           
# disp
disp.summary.df <- data.frame(store = factor(1:88), 
lo = disp.summary[, 1],
md = disp.summary[, 2],
hi = disp.summary[, 3])

disp.summary.df <- disp.summary.df[order(disp.summary.df$md), ]
disp.summary.df$store <- factor(disp.summary.df$store, 
	levels = disp.summary.df$store)

dispCI <- ggplot(disp.summary.df, aes(ymin = lo, ymax = hi, x = store)) + 
geom_linerange() + 
xlab("Store number") +
ylab("Store-level disp term") +
ggtitle("95% Credible intervals for store-level disp term") +
geom_hline(yintercept = 0, linetype = "dashed") +
geom_point(aes(y = md), col = "firebrick3") +
theme_bw() + 
theme(plot.title = element_text(hjust = 0.5))

# logprice.disp
logprice.disp.summary.df <- data.frame(store = factor(1:88), 
lo = logprice.disp.summary[, 1],
md = logprice.disp.summary[, 2],
hi = logprice.disp.summary[, 3])

logprice.disp.summary.df <- logprice.disp.summary.df[order(logprice.disp.summary.df$md), ]
logprice.disp.summary.df$store <- factor(logprice.disp.summary.df$store, 
	levels = logprice.disp.summary.df$store)

logprice.dispCI <- ggplot(logprice.disp.summary.df, aes(ymin = lo, ymax = hi, x = store)) + 
geom_linerange() + 
xlab("Store number") +
ylab("Store-level logprice.disp term") +
ggtitle("95% Credible intervals for store-level logprice.disp term") +
geom_hline(yintercept = 0, linetype = "dashed") +
geom_point(aes(y = md), col = "firebrick3") +
theme_bw() + 
theme(plot.title = element_text(hjust = 0.5))

# Put them all together
mygrid <- grid.arrange(intCI, logpriceCI, dispCI, logprice.dispCI, ncol = 2)

pdf("Cheese/img/CIgrid.pdf", width = 11, height = 7)
grid.arrange(intCI, logpriceCI, dispCI, logprice.dispCI, ncol = 2)
dev.off()

### --------------------------------------------------------------------------
### Credible intervals plot for optimal prices
### --------------------------------------------------------------------------

betamat.nd <- matrix(rep(beta[, 2], each = I), nrow = length.post, byrow = T)
betamat.d <-  matrix(rep(beta[, 4], each = I), nrow = length.post, byrow = T)

# No display
slopes.nd <- b[, seq(2, I*q, by = 4)] + betamat.nd # no display

optprice.nd <- slopes.nd / (slopes.nd + 1)
optprice.nd.summary <- t(as.matrix(apply(optprice.nd, 2, quantile, c(0.025, 0.5, 0.975))))

optprice.nd.summary.df <- data.frame(store = factor(1:88), 
lo = optprice.nd.summary[, 1],
md = optprice.nd.summary[, 2],
hi = optprice.nd.summary[, 3])

optprice.nd.summary.df <- optprice.nd.summary.df[order(optprice.nd.summary.df$md), ]
optprice.nd.summary.df$store <- factor(optprice.nd.summary.df$store, 
	levels = optprice.nd.summary.df$store)

optprice.ndCI <- ggplot(optprice.nd.summary.df, aes(ymin = lo, ymax = hi, x = store)) + 
geom_linerange() + 
xlab("Store number") +
ylab("Store-level optimal price") +
ggtitle("95% Credible intervals for store-level optimal price (no display)") +
geom_point(aes(y = md), col = col0) +
theme_bw() + 
theme(plot.title = element_text(hjust = 0.5))
optprice.ndCI 

# looks bad... remove outliers
remove <- c(which(optprice.nd.summary.df$hi - optprice.nd.summary.df$lo > 5))
optprice.nd.summary.df2 <- optprice.nd.summary.df[-remove, ]
optprice.nd.summary.df2$lo <- optprice.nd.summary.df2$lo * (optprice.nd.summary.df2$lo > 0)

optprice.ndCI2 <- ggplot(optprice.nd.summary.df2, aes(ymin = lo, ymax = hi, x = store)) + 
geom_linerange() + 
xlab("Store number") +
ylab("Store-level optimal price") +
ggtitle("95% Credible intervals for store-level optimal price (no display, outliers removed)") +
geom_point(aes(y = md), col = col0) +
theme_bw() + 
theme(plot.title = element_text(hjust = 0.5))
optprice.ndCI2

### Display

slopes.d <- slopes.nd + b[, seq(4, I*q, by = 4)] + betamat.d

optprice.d <- slopes.d / (slopes.d + 1)
optprice.d.summary <- t(as.matrix(apply(optprice.d, 2, quantile, c(0.025, 0.5, 0.975))))

optprice.d.summary.df <- data.frame(store = factor(1:88), 
lo = optprice.d.summary[, 1],
md = optprice.d.summary[, 2],
hi = optprice.d.summary[, 3])

optprice.d.summary.df <-
optprice.d.summary.df[order(optprice.d.summary.df$md), ]
optprice.d.summary.df$store <- factor(optprice.d.summary.df$store, 
	levels = optprice.d.summary.df$store)

optprice.dCI <- ggplot(optprice.d.summary.df, aes(ymin = lo, ymax = hi, x = store)) + 
geom_linerange() + 
xlab("Store number") +
ylab("Store-level optimal price") +
ggtitle("95% Credible intervals for store-level optimal price (display)") +
geom_point(aes(y = md), col = col1) +
theme_bw() + 
theme(plot.title = element_text(hjust = 0.5))
optprice.dCI 

#### System-wide

cheese.copy <- cheese

optprice.nodisp <- beta[, 2] / (beta[, 2] + 1)
quantile(optprice.nodisp, c(0.025, 0.5, 0.975))

optprice.disp <- (beta[, 2] + beta[4]) / (beta[, 2] + beta[4] + 1)
quantile(optprice.disp, c(0.025, 0.5, 0.975))

optprice.nodisp.hist <- qplot(optprice.nodisp, geom = 'blank') +   
geom_histogram(aes(y = ..density..), 
    alpha = 0.8,  
	bins = 20, 
	colour = "grey95", 
	fill = "grey20") +
xlab("Price / cost ratio") +
ylab("Density") +
ggtitle("Posterior System-wide optimal price / cost ratios (No display)") +
theme_bw() +
theme(legend.position = c(0.85, 0.85),
	plot.title = element_text(hjust = 0.5))
optprice.nodisp.hist 

optprice.disp.hist <- qplot(optprice.disp, geom = 'blank') +   
geom_histogram(aes(y = ..density..), 
    alpha = 0.8,  
	bins = 20, 
	colour = "grey95", 
	fill = "grey20") +
xlab("Price / cost ratio") +
ylab("Density") +
ggtitle("Posterior System-wide optimal price / cost ratios (Display)") +
theme_bw() +
theme(legend.position = c(0.85, 0.85),
	plot.title = element_text(hjust = 0.5))
optprice.disp.hist 



### --------------------------------------------------------------------------
### Create plot of fitted lines
### --------------------------------------------------------------------------

# Posterior mean
b.est <- colMeans(b)
beta.est <- colMeans(beta)

coefest <- matrix(b.est, nrow = I, byrow = T) + 
		   matrix(rep(beta.est, I), nrow = I, byrow = T) 

df.nodisp <- data.frame(store = factor(1:I), 
b0 = coefest[, 1], 
b1 = coefest[, 2])

df.disp <- data.frame(store = factor(1:I), 
b0 = coefest[, 1] + coefest[, 3], 
b1 = coefest[, 2] + coefest[, 4])

# Create plot
allfits <- facetplot + 
geom_abline(data = df.disp, aes(intercept = b0, slope = b1), col = col1) +
geom_abline(data = df.nodisp, aes(intercept = b0, slope = b1), col = "blue") +
ggtitle("Demand curves for all stores with hierarchical model") +
theme(# strip.background = element_blank(),
#        strip.text.x = element_blank(),
	   plot.title = element_text(hjust = 0.5))

pdf("Cheese/img/allfits.pdf", width = 10, height = 10)
allfits
dev.off()

# OLS plot
OLS_HLM <- OLSplot + geom_abline(intercept = beta.est[1], slope = beta.est[2], col = "blue", linetype = 2) + geom_abline(intercept = beta.est[1] + beta.est[3], slope = beta.est[2] + beta.est[4], col = col1, linetype = 2) 
OLS_HLM


# Change order
cheese.copy <- cheese
myorder <- order(df.disp$b1)
df.disp$dispratio <- rep(0, nrow(df.disp))
for (i in 1:nrow(df.disp)) {
	indices <- which(cheese$store == i)
	df.disp$dispratio[i] <- sum(cheese$disp[indices]) / length(indices)
}
myorder <- order(df.disp$dispratio) # order by ratio of data with ads
# myorder <- order(df.disp$b1) # order by slope of display line 
cheese.copy$store <- factor(cheese.copy$store, levels = factor(myorder))

facetplot_order <- ggplot(cheese.copy, aes(logprice, logvol)) + 
geom_point(pch = 1, aes(colour = factor(disp))) + 
xlab("log-Price (USD)") +
ylab("log-Volume sold (units)") +
facet_wrap(~factor(store), ncol = 9) +
scale_colour_manual("Display", values = c("blue", col1)) +
geom_abline(data = df.disp, aes(intercept = b0, slope = b1), col = col1) +
geom_abline(data = df.nodisp, aes(intercept = b0, slope = b1), col = "blue") +
ggtitle("Demand curves for all stores with hierarchical model, ordered") +
theme(plot.title = element_text(hjust = 0.5))

pdf("Cheese/img/facetplot_disporder.pdf", width = 10, height = 10)
facetplot_order
dev.off()
