###########################################################
######### Created by Spencer Woody on 11 Feb 2017 #########
###########################################################

library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer) # display.brewer.all()
library(wesanderson) # nice palettes
library(mvtnorm)

utilities <- read.csv("utilities.csv", header = T)

x <- utilities$temp
y <- log(utilities$gasbill / utilities$billingdays)

source("myfuns03.R")

# ===========================================================================
# Linear smoothers (part one) ===============================================
# ===========================================================================

# nonlinear function f(x)
f1 <- function(x){
	return(x * (x - 4) * (x + 4))
}

# Predictor vector
x1 <- seq(-5, 5, length.out = 40)

# Create sequence along x-space
x.seq <- seq(min(x1), max(x1), length.out = 200)

# Response vector
y1 <- make.noise(x1, f1, "normal", sd = 15)

# Bin width
h1 <- 0.75

# Gaussian kernel smoothing
y.norm <- sapply(
	x.seq, 
	lin.smooth, 
	x = x1, 
	y = y1, 
	kern.fun = kern.norm, 
	h = h1
	)

# Uniform kernel smoothing
y.unif <- sapply(
	x.seq, 
	lin.smooth, 
	x = x1, 
	y = y1, 
	kern.fun = kern.unif, 
	h = h1
	)

# Make a nice plot
h <- qplot(x.seq, geom = "blank") +
xlab("x") +
ylab("y") +
ggtitle(sprintf("Smoothing of cubic function")) +
geom_point(aes(x = x1, y = y1), pch = 1) + 
stat_function(fun = f1, col = col1, linetype = "dashed") + 
geom_line(aes(y = y.norm, colour = "Gaussian kernel")) +
geom_line(aes(y = y.unif, colour = "Uniform kernel")) + 
scale_colour_manual(name = "Smoother", values = c(col3, col2)) +
theme(legend.position = c(0.25, 0.15),
	text = element_text(family="Helvetica"))

h

pdf("firstexample.pdf")
h
dev.off()



# ===========================================================================
# Linear smoothers (cross validation) ==== make big heat map ================
# ===========================================================================

# Sample size
N = 500

# Set limits of x-space
xlo = 0
xhi = 1

# Number of bins for cross-validation
numbins <- 5

# Vector for bandwidths to search over
h.vec <- seq(0.001, 0.125, length.out = 100)

# Vector for standard deviations
s.vec <- seq(0.001, 0.5, length.out = 100)
s.vec <- rev(s.vec)

# Vector for periods
p.vec <- seq(0.1, 1, length.out = 100)

# Matrix for optimal bandwidth values
opt.h <- matrix(nrow = length(p.vec), ncol = length(s.vec)) 

for (i in 1:length(p.vec)) {
	# Select what the period is for this iteration
	p.i <- p.vec[i]
	
	# Create a new function with this current period
	mysin.i <- function(x) {
		return(sin(x * 2*pi / p.i))
	}
	
	for (j in 1:length(s.vec)) {
		# Select what the residual sd is for this iteration
		s.j <- s.vec[j]
		
		# Generate data
		x.ij <- (xhi - xlo) * runif(N) + xlo
		y.ij <- make.noise(x.ij, mysin.i, "normal", sd = s.j)
		
		# Prepare mse matrix for this current iteration
		mse.ij <- matrix(nrow = numbins, ncol = length(h.vec))
		
		# Create random partition into five bins
		jumble <- sample(1:N, N, replace = F)
		bin.indices <- split(jumble, cut(1:N, numbins))
		
		for (bin in 1:numbins) {
			my.indices <- bin.indices[[bin]]
			x.tr <- x.ij[-(my.indices)]
			y.tr <- y.ij[-(my.indices)]
			x.te <- x.ij[my.indices]
			y.te <- y.ij[my.indices]
			mse.ij[bin, ] <- cv(x.tr, y.tr, x.te, y.te, kern.norm, h.vec)
		}
		
		# Average out MSE over all bins
		mse.vec <- colMeans(mse.ij)
		
		# Choose bandwidth with lowest average MSE
		opt.h[i, j] <- h.vec[which.min(mse.vec)]
		
	}
	print(i)
}

# Check to make sure optimal bandwidths look OK
opt.h / max(opt.h)

# Plot these things
ohm <- melt(opt.h)

w <- ggplot(ohm, aes(rev(Var1), Var2, z = value)) +
ggtitle("Picking Optimal Bandwidth for Gaussian Kernel Smoothing (5-fold CV)") + 
xlab("Decreasing standard deviation (0.5:0.001) (less noisy)") +
ylab("Decreasing period (1:0.1) (more wiggly)") +
geom_tile(aes(fill = value)) +
scale_fill_distiller("Bandwidth", palette = "Spectral")
w

# Save this to PDF
pdf("img/opth.pdf", width = 7, height = 6)
w
dev.off()


# ===========================================================================
# Linear smoothers (cross validation) ==== make 2 x 2 plot ==================
# ===========================================================================

# Sample size
N = 500

# Set limits of x-space
xlo = 0
xhi = 1

# Number of bins for cross-validation
numbins <- 5

# Vector for bandwidths to search over
h.vec <- seq(0.001, 0.125, length.out = 100)

# Vector for standard deviations
# s.vec <- seq(0.001, 0.5, length.out = 64)
s.vec <- c(0.1, 0.5)

# Vector for periods
# p.vec <- seq(0.1, 1, length.out = 64)
p.vec <- c(0.125, 1)

# Matrix for optimal bandwidth values
opt.h <- matrix(nrow = length(p.vec), ncol = length(s.vec)) 

# Matrix for random x-values and y-values
x.mat <- matrix(nrow = nrow(opt.h) * ncol(opt.h), ncol = N)
y.mat <- matrix(nrow = nrow(opt.h) * ncol(opt.h), ncol = N)

# Matrix for values of f at x
x.seq <- seq(xlo, xhi, length.out = 200)
fx.mat <- matrix(nrow = nrow(opt.h) * ncol(opt.h), ncol = length(x.seq))
smooth <- matrix(nrow = nrow(opt.h) * ncol(opt.h), ncol = length(x.seq))

count <- 1

for (i in 1:length(p.vec)) {
	# Select what the period is for this iteration
	p.i <- p.vec[i]
	
	# Create a new function with this current period
	mysin.i <- function(x) {
		return(sin(x * 2*pi / p.i))
	}
	
	for (j in 1:length(s.vec)) {
		# Select what the residual sd is for this iteration
		s.j <- s.vec[j]
		
		x.ij <- (xhi - xlo) * runif(N) + xlo
		y.ij <- make.noise(x.ij, mysin.i, "normal", sd = s.j)
		
		fx.mat[count, ] <- mysin.i(x.seq)
		
		x.mat[count, ] <- x.ij
		y.mat[count, ] <- y.ij
		
		mse.ij <- matrix(nrow = numbins, ncol = length(h.vec))
		
		jumble <- sample(1:N, N, replace = F)
		bin.indices <- split(jumble, cut(1:N, numbins))
		
		for (bin in 1:numbins) {
			my.indices <- bin.indices[[bin]]
			x.tr <- x.ij[-(my.indices)]
			y.tr <- y.ij[-(my.indices)]
			x.te <- x.ij[my.indices]
			y.te <- y.ij[my.indices]
			mse.ij[bin, ] <- cv(x.tr, y.tr, x.te, y.te, kern.norm, h.vec)
		}
		
		mse.vec <- colMeans(mse.ij)
		
		# Choose first 4/5 of x's and y's to be training, the rest of testing
		# x.tr <- x.ij[1:round(N*4/5)]
		# y.tr <- y.ij[1:round(N*4/5)]
		# x.te <- x.ij[-(1:round(N*4/5))]
		# y.te <- y.ij[-(1:round(N*4/5))]
		
		# mse.vec <- cv(x.tr, y.tr, x.te, y.te, kern.norm, h.vec)
		
		opt.h[i, j] <- h.vec[which.min(mse.vec)]
		
		smooth[count, ] <- y.norm <- sapply(
							x.seq, 
							lin.smooth, 
							x = x.ij, 
							y = y.ij, 
							kern.fun = kern.norm, 
							h = opt.h[i, j]
							)
							
		
		count <- count + 1
	}
	print(i)
}

y.min <- min(y.mat)
y.max <- max(y.mat)

ess <- qplot(x.seq, geom = "blank") +
xlab("x") +
ylab("y") +
ylim(y.min, y.max) + 
ggtitle(sprintf("sd = %5.3f, T = %5.3f, h = %5.5f", s.vec[1], p.vec[1], opt.h[1, 1])) +
geom_point(aes(x = x.mat[1, ], y = y.mat[1, ]), pch = 1) + 
geom_line(aes(y = fx.mat[1, ]), col = "red", linetype = "dashed") +
geom_line(aes(y = smooth[1, ], colour = "Gaussian kernel")) + 
scale_colour_manual(name = "Smoother", values = "blue") +
theme(legend.position = c(0.25, 0.15), text = element_text(family="Helvetica"))

tee <- qplot(x.seq, geom = "blank") +
xlab("x") +
ylab("y") +
ylim(y.min, y.max) + 
ggtitle(sprintf("sd = %5.3f, T = %5.3f, h = %5.5f", s.vec[2], p.vec[1], opt.h[1, 2])) +
geom_point(aes(x = x.mat[2, ], y = y.mat[2, ]), pch = 1) + 
geom_line(aes(y = fx.mat[2, ]), col = "red", linetype = "dashed") +
geom_line(aes(y = smooth[2, ], colour = "Gaussian kernel")) + 
scale_colour_manual(name = "Smoother", values = "blue") +
theme(legend.position = c(0.25, 0.15), text = element_text(family="Helvetica"))

you <- qplot(x.seq, geom = "blank") +
xlab("x") +
ylab("y") +
ylim(y.min, y.max) + 
ggtitle(sprintf("sd = %5.3f, T = %5.3f, h = %5.5f", s.vec[1], p.vec[2], opt.h[2, 1])) +
geom_point(aes(x = x.mat[3, ], y = y.mat[3, ]), pch = 1) + 
geom_line(aes(y = fx.mat[3, ]), col = "red", linetype = "dashed") +
geom_line(aes(y = smooth[3, ], colour = "Gaussian kernel")) + 
scale_colour_manual(name = "Smoother", values = "blue") +
theme(legend.position = c(0.25, 0.15), text = element_text(family="Helvetica"))

vee <- qplot(x.seq, geom = "blank") +
xlab("x") +
ylab("y") +
ylim(y.min, y.max) + 
ggtitle(sprintf("sd = %5.3f, T = %5.3f, h = %5.5f", s.vec[2], p.vec[2], opt.h[2, 2])) +
geom_point(aes(x = x.mat[4, ], y = y.mat[4, ]), pch = 1) + 
geom_line(aes(y = fx.mat[4, ]), col = "red", linetype = "dashed") +
geom_line(aes(y = smooth[4, ], colour = "Gaussian kernel")) + 
scale_colour_manual(name = "Smoother", values = "blue") +
theme(legend.position = c(0.25, 0.15), text = element_text(family="Helvetica"))
	

pdf("img/2x2.pdf")
grid.arrange(tee, ess, vee, you)
dev.off()

# ===========================================================================
# Local polynomial regression ===============================================
# ===========================================================================

# ---------------------------------------------------------------------------
# Read in the data ----------------------------------------------------------
# ---------------------------------------------------------------------------

utilities <- read.csv("utilities.csv", header = T)

x <- utilities$temp
y <- log(utilities$gasbill / utilities$billingdays)


# ---------------------------------------------------------------------------
# Cross validation ----------------------------------------------------------
# ---------------------------------------------------------------------------

h.vec <- seq(1, 20, length.out = 200)
num.h <- length(h.vec)

hatmat.list <- list()
loocv.vec <- rep(NA, num.h)

for (i in 1:num.h) {
	hatmat.i <- loc.pol.hatmat(x, y, D = 1, h = h.vec[i])
	hatmat.list[[i]] <- hatmat.i
	
	loocv.vec[i] <- loocv(y, hatmat.i)
	
	if ((i %% 20) == 0) {
		print(sprintf("Iteration %i out of %i...", i, num.h))
	}
}

plot(h.vec, loocv.vec, type = "l", xlab = "Bandwidth", ylab = "LOOCV")

h.opt <- h.vec[which.min(loocv.vec)]

# ---------------------------------------------------------------------------
# Make a plot ---------------------------------------------------------------
# ---------------------------------------------------------------------------

x.seq <- seq(min(x), max(x), length.out = 200)

y.smooth <- sapply(
	x.seq,
	loc.pol,
	x.vec = x, 
	y.vec = y,
	D = 1,
	h = h.opt
	)


# Fitted y values
Hatmat <- loc.pol.hatmat(x, y, D = 1, h = h.opt)
y.hat <- Hatmat %*% y

# R-squared
r.sq <- 1 - sum((y - y.hat)^2) / sum((y - mean(y))^2)

# Estimated variance
var.est <- sum((y - y.hat)^2) / 
(length(x) - sum(diag(Hatmat)) + sum(diag(crossprod(Hatmat))))


# Variance of fit
big.hatmat <- matrix(nrow = length(x.seq), ncol = length(y))

fit.var <- rep(NA, nrow(big.hatmat))

for (i in 1:nrow(big.hatmat)) {
	big.hatmat[i, ] <- loc.pol(x.seq[i], x, y, D = 1, h.opt, give.mat = TRUE)$hatmat.vec
	
	fit.var[i] <- var.est * sum((big.hatmat[i, ])^2)
}


# (fix this)
y.lo <- y.smooth - 1.96 * sqrt(fit.var)
y.hi <- y.smooth + 1.96 * sqrt(fit.var)

	
resplot <- qplot(x, y - y.hat, geom = "blank") + 
geom_point(pch = 1) +
geom_hline(yintercept = 0) +
xlab(expression(paste("Temperature (",degree,"F)"))) +
ylab("residual") + 
ggtitle("Residual plot for log-transformed data")

pdf("img/resplot2.pdf")
resplot
dev.off()

q <- qplot(x.seq, geom = "blank") +
xlab(expression(paste("Temperature (",degree,"F)")))  +
ylab("log-Daily gas bill (USD)") +
labs(title = "Daily gas bills for single-family homes in Minnesota") +
geom_ribbon(aes(ymin = y.lo, ymax = y.hi), fill = "grey80") +
geom_point(aes(x = x, y = y), pch = 1) + 
geom_line(aes(y = y.smooth, colour = sprintf("h = %5.4f", h.opt)))  + 
scale_colour_manual(name = "Bandwidth", values = "firebrick3") +
theme(plot.title = element_text(hjust = 0.5), 
text = element_text(family = "Helvetica"),
legend.position = c(0.25, 0.15))

pdf("img/tempplot.pdf")
q
dev.off()

# ===========================================================================
# Gaussian process ==========================================================
# ===========================================================================

x.seq <- seq(0, 1, length.out = 100)

b <- 0.1
tau1.sq <- 1
tau2.sq <- 0

myparams <- c(b, tau1.sq, tau2.sq)

xCM52 <- make.covmat(x.seq, C.M52, params = myparams)
xSE <- make.covmat(x.seq, C.SE, params = myparams)

GP1 <- as.vector(rmvnorm(1, mean = rep(0, 100), sigma = xCM52))
GP2 <- as.vector(rmvnorm(1, mean = rep(0, 100), sigma = xCM52))

tau1.sq.vec <- c(0.001, 0.1, 1e0, 1e1, 1e2)

GP.mat <- matrix(nrow = length(tau1.sq.vec), ncol = length(x.seq))

for (i in 1:length(tau1.sq.vec)) {
	b.i <- 1
	tau1.sq.i <- tau1.sq.vec[i]
	tau2.sq.i <- 1e-6
	
	myparams.i <- c(b, tau1.sq, tau2.sq)
	
	xCM52.i <- make.covmat(x.seq, C.M52, params = myparams.i)
	
	GP.mat[i, ] <- as.vector(rmvnorm(1, mean = rep(0, 100), sigma = xCM52.i))
}	

GP1 <- GP.mat[1, ]
GP2 <- GP.mat[2, ]
GP3 <- GP.mat[3, ]
GP4 <- GP.mat[4, ]
GP5 <- GP.mat[5, ]

cols <- brewer.pal(5, "Set1")


h <- qplot(x.seq, geom = "blank") +
xlab("x") +
ylab("y") +
ggtitle(sprintf("Example of Gaussian processes")) +
geom_line(aes(y = GP1, colour = sprintf("tau1.sq = %f", tau1.sq.vec[1]))) +
geom_line(aes(y = GP2, colour = sprintf("tau1.sq = %f", tau1.sq.vec[2]))) +
geom_line(aes(y = GP3, colour = sprintf("tau1.sq = %f", tau1.sq.vec[3]))) + 
geom_line(aes(y = GP4, colour = sprintf("tau1.sq = %f", tau1.sq.vec[4]))) + 
geom_line(aes(y = GP5, colour = sprintf("tau1.sq = %f", tau1.sq.vec[5]))) +
scale_colour_manual(name = "Params", values = cols) +
theme(legend.position = c(0.25, 0.15),
	text = element_text(family="Helvetica"))

	
	qplot(x.seq, as.vector(GP), geom = "blank") + geom_line()
	plot(x.seq, GP)


h <- qplot(x.seq, geom = "blank") +
xlab("x") +
ylab("y") +
ggtitle(sprintf("Example of Gaussian processes")) +
geom_line(aes(y = GP1, colour = "Gaussian kernel")) +
geom_line(aes(y = GP2, colour = "Uniform kernel")) + 
scale_colour_manual(name = "Params", values = c(col3, col2)) +
theme(legend.position = c(0.25, 0.15),
	text = element_text(family="Helvetica"))



# Plot these things
xCM52.m <- melt(xCM52)
xCM52.m[, 1] <- rep(x.seq, length(x.seq))
xCM52.m[, 2] <- rep(x.seq, each = length(x.seq))

w <- ggplot(xCM52.m, aes(Var1, Var2, z = value)) +
ggtitle("Matern(5/2) covariance function") + 
xlab("x[i]") +
ylab("x[j]") +
geom_tile(aes(fill = value)) +
scale_fill_distiller("Cov(x[i], x[j])", palette = "Spectral") +
theme_bw()
w


x <- c()



# ------------------------------------------------------------------------
# Application to utilities data set --------------------------------------
# ------------------------------------------------------------------------

utilities <- read.csv("utilities.csv", header = T)

x <- utilities$temp
y <- log(utilities$gasbill / utilities$billingdays)

N <- length(x)

b <- 50
tau1.sq <- 5
tau2.sq <- 0

myparams <- c(b, tau1.sq, tau2.sq)

# Choose sigma2
sigma2.init <- 1
pred1 <- GP.predict(x, y, x, C.M52, params = myparams, res.var = sigma2.init)
sigma2.est <- sum((y-pred1$post.mean)^2) / (length(x) - 1)
sigma2.est

# Find tuning parameters

res <- 50

b.vec <- rep(seq(50, 140, length.out = res), res)
tau1.vec <- rep(seq(2, 6, length.out = res), each = res)

ll.vec <- rep(NA, length(b.vec))

for (i in 1:length(ll.vec)) {
	b.i <- b.vec[i]
	tau1.sq.i <- tau1.vec[i]
	tau2.sq.i <- 0

	myparams.i <- c(b.i, tau1.sq.i, tau2.sq.i)
	
	C <- make.covmat(x, C.M52, params = myparams.i)

	ll.vec[i] <-  dmvnorm(y, mean = rep(0, N), 
						  sigma = C + sigma2.est *diag(N), log = TRUE) 
}

tune <- data.frame(cbind(b.vec, tau1.vec, ll.vec))

pdf("img/tune1.pdf", width = 7, height = 6)
ggplot(tune, aes(b.vec, tau1.vec, ll.vec)) + 
geom_tile(aes(fill = ll.vec)) +
geom_contour(aes(z = ll.vec), colour = "white") +
scale_fill_distiller(palette = "Spectral") +
ggtitle("Find tuning parameters for GP on utilities data (log-likelihood)") +
xlab("b") +
ylab("tau1.sq")
dev.off()

which.max(ll.vec)
b.vec[which.max(ll.vec)]
tau1.vec[which.max(ll.vec)]

opt.params <- c(b.vec[which.max(ll.vec)], tau1.vec[which.max(ll.vec)], 0)

# Now make final fit 

x.seq <- seq(min(x), max(x), length.out = 200)

pred2 <- GP.predict(x, y, x.seq, C.M52, params = myparams, res.var = sigma2.est)

fit <- pred2$post.mean

f.lo <- fit - 1.96 * sqrt(pred2$post.var)
f.hi <- fit + 1.96 * sqrt(pred2$post.var)

GP1 <- qplot(x.seq, geom = "blank") +
xlab(expression(paste("Temperature (",degree,"F)")))  +
ylab("log-Daily gas bill (USD)") +
labs(title = "Daily gas bills for single-family homes in Minnesota") +
geom_ribbon(aes(ymin = f.lo, ymax = f.hi), fill = "grey80") +
geom_point(aes(x = x, y = y), pch = 1) + 
geom_line(aes(y = fit,colour = sprintf("b=%.3f, tau1.sq=%.3f", opt.params[1], opt.params[2])))  + 
theme_bw() +
scale_colour_manual(name = "GP with M52", values = "firebrick3") +
theme(plot.title = element_text(hjust = 0.5), 
text = element_text(family = "Helvetica"),
legend.position = c(0.25, 0.15)) 
GP1

pdf("img/Minn_GP.pdf")
GP1
dev.off()

# Find optimal parameters

# ===========================================================================
# Weather data ==============================================================
# ===========================================================================



weather <- read.csv("weather.csv", header = T)

# Prep the color palette
cols <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100))

# Temperature

# Plot the data
pplot <- ggplot(weather, aes(x=lon, y=lat, colour = pressure)) + 
geom_point(pch = 15, size = 2) + 
scale_colour_gradientn(colours=cols) +
theme_bw() 

pdf("img/pplot.pdf", width = 7, height = 5)
pplot
dev.off()

tplot <- ggplot(weather, aes(x=lon, y=lat, colour = temperature)) + 
geom_point(pch = 15, size = 2) + 
scale_colour_gradientn(colours=cols) +
theme_bw() 

pdf("img/tplot.pdf", width = 7, height = 5)
tplot
dev.off()

# Choose sigma2 and hyperparameters
myparams <- c(20, 5, 0)
dist.params <- c(sqrt(2), 1)

sigma2.init <- 3
wpred1 <- GP.predict2(weather[, 3:4], weather[, 2], weather[, 3:4], C.M52, params = myparams, dist.params = c(1,2), res.var = sigma2.init)

sigma2.est2 <- sum((weather[, 2] - wpred1$post.mean)^2) / (nrow(weather) - 1)
sigma2.est2

# Choose tuning parameters

N1 <- 30

b.vec <- rep(seq(0.5, 3, length.out = N1), N1)
tau1.vec <- rep(seq(3, 7, length.out = N1), each = N1)

ll.vec <- rep(NA, length(b.vec))


for (i in 1:length(ll.vec)) {
	b.i <- b.vec[i]
	tau1.sq.i <- tau1.vec[i]
	tau2.sq.i <- 0

	myparams.i <- c(b.i, tau1.sq.i, tau2.sq.i)
	
	C <- make.covmat2(weather[, 3:4], C.M52, params = myparams.i, dist.params = dist.params <- c(1 / sqrt(2), 1))

	ll.vec[i] <-  dmvnorm(weather[, 2], mean = rep(0, nrow(weather)), 
					sigma = C + sigma2.est2 * diag(nrow(weather)), log = TRUE) 
}

tune2 <- data.frame(cbind(b.vec, tau1.vec, ll.vec))

ggplot(tune2, aes(b.vec, tau1.vec, ll.vec)) + 
geom_tile(aes(fill = ll.vec)) +
geom_contour(aes(z = ll.vec), colour = "white") +
scale_fill_distiller(palette = "Spectral") 

which.max(ll.vec)

# Make plot of smooth version


# mapImageData1 <- get_map(location = c(lon = -0.016179, lat = 51.538525),
#     color = "color",
#     source = "google",
#     maptype = "satellite",
#     zoom = 17)
#
# ggmap(mapImageData1,
#     extent = "device",
#     ylab = "Latitude",
#     xlab = "Longitude")
#
# 	map <- get_map("Atlanta", zoom=10)
# 	p <- ggmap(map)


myparams <- c(b.vec[which.max(ll.vec)], tau1.vec[which.max(ll.vec)], 0)


NN <- 64

lon.vec <- rep(seq(min(weather$lon), max(weather$lon), length.out = NN), NN)
lat.vec <- rep(seq(min(weather$lat), max(weather$lat), length.out = NN), each = NN)

new <- data.frame(cbind(lon.vec, lat.vec))
names(new)
names(new) <- names(weather)[c(3, 4)]

fit.vec <- GP.predict2(weather[, 3:4], weather[, 2], new, C.M52, params = myparams, dist.params = c(1,2), res.var = sigma2.est2)$post.mean


wfit <- data.frame(cbind(lon.vec, lat.vec, fit.vec))

pdf("img/GPtemp.pdf", width = 8, height = 6)
ggplot(wfit, aes(lon.vec, lat.vec, fit.vec)) + 
geom_tile(aes(fill = fit.vec)) + 
geom_contour(aes(z = fit.vec), colour = "white") +
scale_fill_distiller(name = "temp", palette = "RdYlBu") +
theme_bw() +
ggtitle("Fitted GP for temperature data")
dev.off()




cols <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100))

bounds <- c(min(weather$lon), min(weather$lat), max(weather$lon), max(weather$lat))

bw.map <- get_map(location = bounds, maptype = "satellite")

ggmap(bw.map) +
geom_point(data = weather, aes(x = lon, y = lat, colour = temperature), size = 2, pch = 19) + 
scale_fill_gradientn(colours=cols) 

ggmap(bw.map) + geom_point(data = weather, aes(x = lon, y = lat, fill = temperature), colour = 'gray', size = 2, pch = 21) + 
scale_fill_gradientn(colours=cols) 

ggmap(bw.map) +
  geom_tile(data = pred.values, aes(x = lon, y = lat, fill = pred), alpha = 0.6) +
  scale_fill_distiller(palette = 'RdYlBu') +
  geom_contour(data = pred.values, aes(x = lon, y = lat, z = pred, colour = ..level..), size = 1) +
  scale_colour_gradientn(colors = cols) +
  geom_point(data = weather, aes(x = lon, y = lat, fill = pressure), colour = 'gray', size = 2, pch = 21)

ggmap(bw.map) +
  geom_contour(data = pred.values, aes(x = lon, y = lat, z = pred, colour = ..level..), size = 1) +
  scale_colour_gradientn(colors = cols)






# ===========================================================================
# Extra code for CV plotting.... ============================================
# ===========================================================================

# Cross-validation with error bars

RSS.vec <- apply(RSS.mat, 2, mean)
RSS.SE <- apply(RSS.mat, 2, sd) / sqrt(numbins)

lower = RSS.vec - RSS.SE
upper = RSS.vec + RSS.SE

# Make a plot!
pdf("perror.pdf", width = 12 / 1.25, height = 8 / 1.25)
h <- qplot(log(fit1$lambda), RSS.vec, geom = "path") 
h + xlab(expression(paste("Penalization term, log(", lambda, ")"))) + 
ylab("Expected prediction error") + 
labs(title = sprintf("%i-fold Cross-validation, Mallow's CP, and In-sample MSE", numbins)) +
geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80") + 
geom_line(aes(y = RSS.vec, colour = "CV"), show.legend = TRUE) +
geom_line(aes(y = MSE.lasso, colour = "In-sample MSE"), show.legend = TRUE) +
geom_line(aes(y = MCP, colour = "CP"), show.legend = TRUE) + 
geom_vline(xintercept = log(minlambdaCV), linetype = 3, col = "black", size = 0.75, show.legend = TRUE) +
geom_vline(xintercept = log(minlambdaMCP), linetype = 3, col = "red", size = 0.75, show.legend = TRUE) +
scale_colour_manual(name = "", values = c("CV" = "black", "In-sample MSE" = "blue", "CP" = "red")) 
dev.off()


numbins <- 10
jumble <- sample(1:N2, N2, replace = F)
bin.indices <- split(jumble, cut(1:N2, numbins))











