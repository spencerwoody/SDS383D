###########################################################
######### Created by Spencer Woody on 27 Mar 2017 #########
###########################################################

library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer) # display.brewer.all()
library(xtable)
library(lme4)
library(mosaic)
library(dplyr)

# Read in the data

cheese.csv <- read.csv("Cheese/cheese.csv")

cheeseplot <- ggplot(cheese.csv, aes(log(price), log(vol))) + 
geom_point(pch = 1, aes(color = factor(disp))) +
scale_colour_manual("Display", values = c("blue", "orange")) 


cheese <- data.frame(store = factor(cheese.csv$store),
	logprice = log(cheese.csv$price),
	disp = cheese.csv$disp,
	logprice.disp = log(cheese.csv$price) * cheese.csv$disp,
	logvol = log(cheese.csv$vol))

stores.list <- unique(cheese$store)

levels(cheese$store) <- 1:length(levels(cheese$store))
cheese$store <- as.numeric(cheese$store)

cheese$disp <- as.numeric(cheese$disp)

# Demand curve for all stores (individually)
# Small multiples
# ... not very helpful
myplot <- ggplot(cheese, aes(x, y)) + 
geom_point(pch = 1, aes(colour = z)) + 
geom_abline(slope = betahat[2] + betahat[4], intercept = betahat[1] + betahat[3]) +
xlab("log-Price (USD)") +
ylab("log-Volume sold (units)") +
ggtitle("Demand curve for all stores") +
facet_wrap(~store, ncol = 9) +
scale_colour_manual("Display", values = c("blue", "orange")) +
theme(strip.background = element_blank(),
       strip.text.x = element_blank(),
	   plot.title = element_text(hjust = 0.5))




myplot <- ggplot(cheese, aes(logprice, logvol)) + 
geom_point(pch = 1, aes(colour = factor(disp))) + 
xlab("log-Price (USD)") +
ylab("log-Volume sold (units)") +
ggtitle("Demand curve for all stores") +
facet_wrap(~factor(store), ncol = 9) +
scale_colour_manual("Display", values = c("blue", "orange")) +
theme(strip.background = element_blank(),
       strip.text.x = element_blank(),
	   plot.title = element_text(hjust = 0.5))









pdf("Cheese/img/panel.pdf", height = 10, width = 10)
p
dev.off()

# Demand curve for all stores (pooled)
q <- ggplot(cheese, aes(logvol, logprice)) + 
geom_abline(slope = betahat[2] + betahat[4], intercept = betahat[1] + betahat[3]) + 
xlab("log-Price (USD)") +
ylab("log-Volume sold (units)") +
ggtitle("Demand curve for all stores") +
geom_point(pch = 1, aes(colour = factor(disp))) +
scale_colour_manual("Display", values = c("blue", "orange")) +
# theme_bw() +
theme(plot.title = element_text(hjust = 0.5)) 

pdf("Cheese/img/cheesescatter.pdf")
q
dev.off()

# Fit an HLM
my.hlm <- lmer(y ~ x + z + z * x + (1 + x + z + z * x | store), data = cheese)
summary(my.hlm)

my.hlm <- lmer(y ~ x + z + z * x + (1 + x + z + z * x | store), data = cheese) 

# Random effects
ran <- ranef(my.hlm, condVar = T)
dotplot(ran)
ggCaterpillar(ran)

# Residual plot
resid.plot <- qplot(fitted(my.hlm), resid(my.hlm), geom = "blank") + 
geom_hline(yintercept = 0) +
geom_point(pch = 1, col = "dodgerblue3") +
ggtitle("Residual plot for HLM of Cheese Demand") +
xlab("Fitted value of log-Volume sold (units)") +
ylab("Residual") +
# theme_bw() +
theme(plot.title = element_text(hjust = 0.5)) 

pdf("Cheese/img/cheeseresid1.pdf")
resid.plot
dev.off()

# http://stackoverflow.com/questions/13847936/in-r-plotting-random-effects-from-lmer-lme4-package-using-qqmath-or-dotplot

ggCaterpillar <- function(re, QQ=TRUE, likeDotplot=TRUE) {
    require(ggplot2)
    f <- function(x) {
        pv   <- attr(x, "postVar")
        cols <- 1:(dim(pv)[1])
        se   <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
        ord  <- unlist(lapply(x, order)) + rep((0:(ncol(x) - 1)) * nrow(x), each=nrow(x))
        pDf  <- data.frame(y=unlist(x)[ord],
                           ci=1.96*se[ord],
                           nQQ=rep(qnorm(ppoints(nrow(x))), ncol(x)),
                           ID=factor(rep(rownames(x), ncol(x))[ord], levels=rownames(x)[ord]),
                           ind=gl(ncol(x), nrow(x), labels=names(x)))

        if(QQ) {  ## normal QQ-plot
            p <- ggplot(pDf, aes(nQQ, y))
            p <- p + facet_wrap(~ ind, scales="free")
            p <- p + xlab("Standard normal quantiles") + ylab("Random effect quantiles")
        } else {  ## caterpillar dotplot
            p <- ggplot(pDf, aes(ID, y)) + coord_flip()
            if(likeDotplot) {  ## imitate dotplot() -> same scales for random effects
                p <- p + facet_wrap(~ ind)
            } else {           ## different scales for random effects
                p <- p + facet_grid(ind ~ ., scales="free_y")
            }
            p <- p + xlab("Levels") + ylab("Random effects")
        }

        p <- p + theme(legend.position="none")
        p <- p + geom_hline(yintercept=0)
        p <- p + geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0, colour="black")
        p <- p + geom_point(aes(size=1.2), colour="blue") 
        return(p)
    }

    lapply(re, f)
}

#####



library(ggplot2)
library(reshape2)
library(gridExtra)
library(mvtnorm)
library(sparseMVN)
library(RColorBrewer) # display.brewer.all()
library(xtable)
library(lme4)
library(mosaic)
library(dplyr)

# Read in the data

cheese.csv <- read.csv("Cheese/cheese.csv")

cheese <- data.frame(store = factor(cheese.csv$store),
	logprice = log(cheese.csv$price),
	disp = cheese.csv$disp,
	logprice.disp = log(cheese.csv$price) * cheese.csv$disp,
	logvol = log(cheese.csv$vol))

stores.list <- unique(cheese$store)

levels(cheese$store) <- 1:length(levels(cheese$store))
cheese$store <- as.numeric(cheese$store)





niter <- 5000
nburn <- 1000





N <- nrow(cheese)

cheese.mat <- as.matrix(cheese)

inds <- cheese.mat[, 1]
X <- cbind(rep(1, N), cheese.mat[, 2:4])
Z <- cbind(rep(1, N), cheese.mat[, 2:4])
y <- cheese.mat[, 5]

I <- length(unique(inds))

p <- ncol(X)
q <- ncol(X)

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


SSX <- matrix(0, nrow = p, ncol = p) 

for (i in 1:I) {
	SSX <- SSX + XtX.block[((i - 1) * p + 1):(i * p), ((i - 1) * p + 1):(i * p)]
	# SSX <- SSX + crossprod(X[inds == i, ])
}

SSX.inv <- matrix(solve(SSX), ncol = p)

b <- matrix(nrow = niter, ncol = I * q)
beta <- matrix(nrow = niter, ncol = p)
lambda <- rep(NA, niter)
D <- array(dim = c(p, p, niter))

nu <- 5
Psi <- diag(q)

b[1, ] <- rep(0.1, q * I)
beta[1, ] <- rep(0.1, p)

# b[1, ] <- rep((as.numeric(betahat) / 2), I)
# beta[1, ] <- as.numeric(betahat) / 2

D[, , 1] <- diag(q)
lambda[1] <- 1


for (iter in 2:niter) {
	b.c <- b[iter - 1, ]
	lambda.c <- lambda[iter - 1]
	beta.c <- beta[iter - 1, ]
	D.c <- D[, , iter - 1]
	
	#########
	# beta.c <- rep(0, p)
	#########
	
	
	# Update lambda
	RSS <- sum((y - X %*% beta.c - Z.block %*% b.c)^2)
	lambda.n <- rgamma(1, N / 2, RSS / 2)
	
	# Update b
	v <- y - X %*% beta.c
	
	b.prec <- lambda.n * ZtZ.block + bdiag(rep(list(D.c), I))
	b.mean <- lambda.n * solve(b.prec, crossprod(Z.block, v))

	mychol <- Cholesky(b.prec)

	b.n <- as.numeric(rmvn.sparse(1, b.mean, mychol, TRUE))
	
	# b.cov <- solve(lambda.n * ZtZ.block + bdiag(rep(list(D.c), I)) )
	# b.mean <- lambda.n * b.cov %*% crossprod(Z.block, v)
	#
	# b.n <- as.numeric(rmvnorm(1, b.mean, matrix(b.cov, nrow = I * q)))
	
	# Update D
	B.n <- matrix(as.numeric(b.n), nrow = q)
	D.n <- rWishart(1, nu + I, Psi + tcrossprod(B.n))
	
	# Update beta
	w <- y - Z.block %*% b.n
	beta.n <- rmvnorm(1, SSX.inv %*% crossprod(X, w), SSX.inv / lambda.n)
	
	if (iter %% 100 == 0) {
		print(sprintf("Iteration %i out of %i...", iter, niter))
	}
	
	# print(iter)
	
	# Add updates to vectors / matrices / arrays
	b[iter, ] <- b.n
	lambda[iter] <- lambda.n
	beta[iter, ] <- beta.n
	D[, , iter] <- D.n
}

b <- b[-(1:nburn), ]
lambda <- lambda[-(1:nburn)]
beta <- beta[-(1:nburn), ]
D <- D[, , -(1:nburn)]

plot(b[, 1], type = 'l')

plot(lambda, type = "l")


b.est <- colMeans(b)
beta.est <- colMeans(beta)

matrix(rep(beta.est, I), nrow = I, byrow = T) 

coefest <- matrix(b.est, nrow = I, byrow = T) + matrix(rep(beta.est, I), nrow = I, byrow = T) 

df1 <- data.frame(store = factor(unique(inds)), b0 = coefest[, 1], b1 = coefest[, 2])
df2 <- data.frame(store = factor(unique(inds)), 
b0 = coefest[, 1] + coefest[, 3], 
b1 = coefest[, 2] + coefest[, 4])


df1 <- data.frame(store = factor(unique(inds)), b0 = betahat[1], b1 = betahat[2])
df2 <- data.frame(store = factor(unique(inds)), b0 = betahat[1] + betahat[3], b1 = betahat[2] + betahat[4])



myplot + 
geom_abline(data = df1, aes(intercept = b0, slope = b1), col = "blue") +
geom_abline(data = df2, aes(intercept = b0, slope = b1), col = "orange")

myplot + geom_abline(intercept = betahat[1], slope = betahat[2], col = "blue")

beta.c <- matrix(rep(0.1, 4), ncol = 1)
b.c <- matrix(rep(0.1, 4 * I), ncol = 1)



# Gibbs sampler

hlm.gibbs <- function(y, X, Z, inds, init.par, hyper) {
	X.list <- 
	Z.list <- list()
	
}
