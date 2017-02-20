###########################################################
######### Created by Spencer Woody on 11 Feb 2017 #########
###########################################################

library(ggplot2)
library(wesanderson) # nice palettes

# Prep color palette
pal <- wes_palette("Zissou", 5)
col1 <- pal[5]
col2 <- pal[4]
col3 <- pal[1]

col1 <- "red"
col2 <- "orange"
col3 <- "blue"

source("myfuns03.R")

# ===========================================================================
# Linear smoothers ==========================================================
# ===========================================================================

period <- 0.25

f1 <- function(x){
	return(sin(x * 2*pi / period ))
}

# # nonlinear function f(x)
# f1 <- function(x){
# 	return(x * (x - 4) * (x + 4))
# }

# Predictor vector
x1 <- seq(-5, 5, length.out = 40)

# Sample size
N = 200

# Set limits of x-space
xlo = 0
xhi = 1

# Draw uniformly across this space
x1 <- (xhi - xlo) * runif(500) + xlo

# Create sequence along x-space
x.seq <- seq(xlo, xhi, length.out = 200)

# Response vector
y1 <- make.noise(x1, f1, "normal", sd = 1/2)

# Bin width
h1 <- 0.01436364

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
xlab(expression(paste("Temperature (",degree,"C)"))) +
ylab("y") +
ggtitle("Smoothing functions") +
geom_point(aes(x = x1, y = y1), pch = 1) + 
stat_function(fun = f1, col = col1, linetype = "dashed") + 
geom_line(aes(y = y.norm, colour = "Gaussian kernel")) +
geom_line(aes(y = y.unif, colour = "Uniform kernel")) + 
scale_colour_manual(name = "Smoother", values = c(col3, col2)) +
theme(legend.position = c(0.15, 0.15),
	text = element_text(family="Helvetica"))

h

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


jumble <- sample(1:length(x1), round(length(x1) / 5), rep = F)

x.tr <- x1[-jumble]
y.tr <- y1[-jumble]
x.te <- x1[jumble]
y.te <- y1[jumble]

h1 <- seq()

y.pr <- sapply(
	x.te, 
	lin.smooth, 
	x = x.tr, 
	y = y.tr, 
	kern.fun = kern.unif, 
	h = h1
	)

mean((y.pr - y.te)^2)

h.vec <- seq(0.001, 0.05, length.out = 100)

mse.vec <- cv(x.tr, y.tr, x.te, y.te, kern.norm, h.vec)

which.min(mse.vec)
h.vec[which.min(mse.vec)]

plot(h.vec, mse.vec, type = "l")

# ===========================================================================
# Gaussian process ==========================================================
# ===========================================================================

x.seq <- seq(0, 1, length.out = 100)

b <- 1
tau1.sq <- 1e-6
tau2.sq <- 1e-5

myparams <- c(b, tau1.sq, tau2.sq)

xCM52 <- make.covmat(x.seq, C.M52, params = myparams)
xSE <- make.covmat(x.seq, C.SE, params = myparams)














