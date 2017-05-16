### ==========================================================================
### --------------------------------------------------------------------------
### Created by Spencer Woody on 06 May 2017
### --------------------------------------------------------------------------
### ==========================================================================

### --------------------------------------------------------------------------
### Libraries / source code
### --------------------------------------------------------------------------

library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(Matrix)
library(mvtnorm)
library(BayesLogit) # Sample from Polya-Gamma
library(DESeq2) # Normalize data, give EB estimates of overdispersion parameters

# If DESeq2 is not installed, use these two commands
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

# Covariance matrix functions for GP
source("covfuns.R")

### --------------------------------------------------------------------------
### Data prep and dispersion estimates
### --------------------------------------------------------------------------

# Vector of times when samples are collected
times <- c(3, 4, 5, 6, 8, 24, 48, 168, 336)

# Read in raw data 
# Be sure to remove rows with all zeros later...
RNA <- read.csv("ecoli_RNA.csv", header = T	)

# As a matrix
RNAmat <- as.matrix(RNA[, -1])

# Obtain EB estimates of overdispersions and size factors from DESeq2
replicate <- factor(c(rep("rep1", 9), rep("rep2", 9), rep("rep3", 9)))
coldata <- data.frame(row.names=colnames(RNA)[-1], replicate)

dds <- DESeqDataSetFromMatrix(
	countData=RNAmat, 
	colData=coldata, 
	design = ~replicate)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

a   <- 1 / dispersions(dds) # reciprocal is important: see notation
ssf <- sizeFactors(dds)

# Normalize the matrix by size factors
RNAnorm <- round(sweep(RNAmat, 2, ssf, `/`))

# Example plots
myind <- 10

par(mfrow = c(1, 3))
plot(RNAnorm[myind, 1:9], type = "l", ylim = c(0, max(RNAnorm[myind, ])))
points(RNAnorm[myind, 1:9])

plot(RNAnorm[myind, 10:18], type = "l", ylim = c(0, max(RNAnorm[myind, ])))
points(RNAnorm[myind, 10:18])

plot(RNAnorm[myind, 19:27], type = "l", ylim = c(0, max(RNAnorm[myind, ])))
points(RNAnorm[myind, 19:27])


### --------------------------------------------------------------------------
### Gibbs sampler [[ Everything but MH step for now ... ]]
### --------------------------------------------------------------------------

# Choose a gene
mygene <- 2214

# Data vector
y.n <- as.numeric(RNAnorm[mygene, ])
a.n <- a[mygene]

z.n <- (y.n - a.n) / 2

# Time vector
t.n <- log(rep(times, 3))
t.nr <- log(times)
NN <- length(t.n) # Total number of samples

# Replicate
reps <- rep(c("rep1", "rep2", "rep3"), each = length(times))
numreps <- length(unique(reps))

# Create a dataframe
gene.data <- data.frame(replicate = reps, time = t.n, read = y.n)

ggplot(gene.data, aes(x = time, y = read)) + 
geom_point(aes(col = replicate)) +
geom_line(aes(col = replicate))

# Parameters of covariance functions
theta.g   <- c(2, 1, 0)
theta.p   <- c(2, 1, 0)


## STEP 0 - Update initialize latents, create cov matrix ---------------------

# (sample from prior)
omega.n <- rpg(num = NN, h = y.n + a.n, z = rep(0, NN))
Omega.n <- Diagonal(x = omega.n)

## Create covariance matrix, precache its inverse ----------------------------
K.g <- make.covmat(t.n, cov.fun = C.M52, params = theta.g)
K.p <- make.covmat(t.nr, cov.fun = C.M52, params = theta.p)

K.n <- K.g + bdiag(rep(list(K.p), 3))
K.nINV <- solve(K.n)

## STEP 1 - Update psi -------------------------------------------------------

Sigma.n <- as.matrix(solve(K.nINV + Omega.n))
psi.n <- as.numeric(rmvnorm(1, Sigma.n %*% z.n, Sigma.n))

## STEP 2 - Sample latent variables ------------------------------------------

omega.n <- rpg(num = NN, h = y.n + a.n, z = psi.n)
Omega.n <- Diagonal(x = omega.n)

## STEP 3 - Sample prediction ------------------------------------------------

# How many points to smooth over?
pred.length <- 100

# Time vector to sample over
t.n.star <- seq(min(t.n), max(t.n), length.out = pred.length)






## STEP 4 - Find start for params --------------------------------------------
#
# par.init <- c(5, 1, 2, 1)
#
# myoptim <- optim(par = log(par.init), fn = neg.ll, psi.n = psi.n, z.n = z.n)
#
# myoptim
# PAR <- exp(myoptim$par)
# PAR
#
# theta.g <- c(PAR[1], PAR[2], 0)
# theta.p <- c(PAR[3], PAR[3], 0)
#



## Create objects for Gibbs sampler ------------------------------------------

niter <- 4e3
nburn <- 1e3

omega.n_keep <- matrix(nrow = niter, ncol = NN)
psi.n_keep <- matrix(nrow = niter, ncol = NN)
g.n.star_keep <- matrix(nrow = niter, ncol = pred.length)
psi.nr.star_keep <- array(dim = c(pred.length, numreps, niter))

for (iter in 1:(nburn + niter)) {
	
	### ------------------------------------------------------------
	### Update psi  ------------------------------------------------
	### ------------------------------------------------------------
	
	Sigma.n <- as.matrix(solve(K.nINV + Omega.n))
	psi.n <- as.numeric(rmvnorm(1, Sigma.n %*% z.n, Sigma.n))
	
	### ------------------------------------------------------------
	### Update latent variables  -----------------------------------
	### ------------------------------------------------------------
	
	# Update latent variable
	omega.n <- rpg(num = NN, h = y.n + a.n, z = psi.n)
	Omega.n <- Diagonal(x = omega.n)
	
	### ------------------------------------------------------------
	### Step for sampling hyperparameters... -----------------------
	### ------------------------------------------------------------
	
	# Update theta.g and theta.p ...
	
	#
	# MH step
	
	# x.p <- x * exp(tune1 * rnorm(1)) 
	# alpha <- likep / like * x.p / x
	
	#
	
	# Update K.nINV
	
	# K.g <- make.covmat(t.n, cov.fun = C.M52, params = theta.g)
	# K.p <- make.covmat(t.nr, cov.fun = C.M52, params = theta.p)
	#
	# K.n <- K.g + bdiag(rep(list(K.p), 3))
	# K.nINV <- solve(K.n)
	
	### ------------------------------------------------------------
	### Prediction step --------------------------------------------
	### ------------------------------------------------------------
	
	# Gene-level function
	K.n.star <- make.covmat(t.n.star, t.n, cov.fun = C.M52, params = theta.g)

	K.n.starstar <- make.covmat(t.n.star, cov.fun = C.M52, params = theta.g)

	g.n.star.mean <- K.n.star %*% K.nINV %*% Sigma.n %*% z.n
	g.n.star.cov  <- as.matrix(
		K.n.star %*% K.nINV %*% Sigma.n %*% K.nINV %*% t(K.n.star) +
		K.n.starstar - K.n.star %*% K.nINV %*% t(K.n.star)
		)

	g.n.star <- as.numeric(rmvnorm(1, g.n.star.mean, g.n.star.cov))
	
	# Replicate-level functions
	psi.nr.star <- matrix(nrow = pred.length, ncol = numreps)
	
	K.nr.starstar <- K.n.starstar + 
	make.covmat(t.n.star, cov.fun = C.M52, params = theta.p)
	
	smallmat <- make.covmat(t.n.star, t.nr, cov.fun = C.M52, params = theta.p)
	
	for (r in 1:numreps) {
		indices <- ((r - 1) * 9 + 1):(r * 9)
		
		K.nr.star <- K.n.star
		K.nr.star[, indices] <- K.nr.star[, indices] + smallmat
		
		psi.nr.star.mean <- K.nr.star %*% K.nINV %*% Sigma.n %*% z.n
		psi.nr.star.cov  <- as.matrix(
			K.nr.star %*% K.nINV %*% Sigma.n %*% K.nINV %*% t(K.nr.star) +
			K.nr.starstar - K.nr.star %*% K.nINV %*% t(K.nr.star)
			)
		
		psi.nr.star[, r] <- as.numeric(rmvnorm(
			1, 
			psi.nr.star.mean, 
			psi.nr.star.cov))
	}
	
	### ------------------------------------------------------------
	### Store update -----------------------------------------------
	### ------------------------------------------------------------
	
	# Store updates
	if (iter > nburn) {
		omega.n_keep[iter - nburn, ] <- omega.n
		psi.n_keep[iter - nburn, ] <- psi.n
		g.n.star_keep[iter - nburn, ] <- g.n.star
		psi.nr.star_keep[, , iter - nburn] <- psi.nr.star
	}
	
	### ------------------------------------------------------------
	### Print message ----------------------------------------------
	### ------------------------------------------------------------
	
	if (iter %% 100 == 0) {
		print(sprintf("Iteration %i out of %i", iter, (niter + nburn)))
	} 
}

g.post <- t(as.matrix(apply(g.n.star_keep, 2, quantile, c(0.025, 0.5, 0.975))))

p.post <- apply(psi.nr.star_keep, c(1, 2), quantile, c(0.025, 0.5, 0.975))

p.post1 <- t(as.matrix(p.post[, , 1]))
p.post2 <- t(as.matrix(p.post[, , 2]))
p.post3 <- t(as.matrix(p.post[, , 3]))

g.df <- data.frame(
	time = t.n.star, 
	lo = g.post[, 1], 
	md = g.post[, 2], 
	hi = g.post[, 3]
	)

p.df <- data.frame(
	replicate = rep(unique(reps), each = pred.length),
	time = rep(t.n.star, 3),
	lo = c(p.post1[, 1], p.post2[, 1], p.post3[, 1]),
	md = c(p.post1[, 2], p.post2[, 2], p.post3[, 2]),
	hi = c(p.post1[, 3], p.post2[, 3], p.post3[, 3])
	)

ylo <- min(c(p.df$lo, g.df$lo))
yhi <- max(c(p.df$hi, g.df$hi))

# plot(omega.n_keep[, 1], type = "l")
# plot(psi.n_keep[, 1], type = "l")

# gene level function
gene_plot <- ggplot(g.df, aes(x= time, y = md)) + 
geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey", alpha = 0.6) + 
geom_line(aes(x = time, y = md)) + 
geom_line(data = p.df, aes(x = time, y = md, col = replicate), lty = "dashed") +
ylab("Time course function") +
xlab("log-Time (hr)") + 
ggtitle("Gene-level") +
ylim(ylo, yhi) +
theme_bw() +
theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
gene_plot

# replicate level functions
rep_plot <- ggplot(p.df, aes(x= time, y = md)) + 
geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey", alpha = 0.6) + 
geom_line(aes(x = time, y = md)) + 
geom_point(data = gene.data, aes(x = time, y = log(read / a.n), col = replicate)) +
facet_wrap(~replicate) +
ylab("Time course function") +
xlab("log-Time (hr)") +
ggtitle(sprintf("Replicate-level (Gene %i)", mygene)) +
ylim(ylo, yhi) +
theme_bw() +
theme(legend.justification=c(0.975,0.025), legend.position=c(0.975,0.025), legend.direction = "horizontal",   strip.text.x = element_blank(),
plot.title = element_text(hjust = 0.5), legend.background = element_rect(fill="gray95", size=1, linetype="dotted")) 
rep_plot

# Save the plot to a PDF
pdf(sprintf("img/grid_gene%i.pdf", mygene), height = 2.5, width = 8.5)
grid.arrange(gene_plot, rep_plot, ncol = 2, widths = c(1.1, 3))
dev.off()

# Plot the data
ggplot(gene.data, aes(x = time, y = read)) + 
geom_point(aes(col = replicate)) +
geom_line(aes(col = replicate))






# exponentiate
gene_plot2 <- ggplot(g.df, aes(x= time, y = a.n * exp(md))) + 
geom_ribbon(aes(ymin = a.n * exp(lo), ymax = a.n * exp(hi)), fill = "grey") + 
geom_line(aes(x = time, y = a.n * exp(md)))

rep_plot2 <- ggplot(p.df, aes(x= time, y = a.n * exp(md))) + 
geom_ribbon(aes(ymin = a.n * exp(lo), ymax = a.n * exp(hi)), fill = "grey") + geom_line(aes(x = time, y = a.n * exp(md), col = replicate)) + 
geom_point(data = gene.data, aes(x = time, y = read, col = replicate)) +
facet_wrap(~replicate)

grid.arrange(gene_plot2, rep_plot2, ncol = 2)





mysample <- sample(1:nrow(RNAnorm), 1, replace = F)
mysample

# RNAnorm[sample(1:nrow(RNAnorm), 4, replace = F), 1:9]


# Choose a gene

mysample <- sample(1:nrow(RNAnorm), 1, replace = F)
mysample

mygene2 <- mysample

# Data vector
y.n2 <- as.numeric(RNAnorm[mygene2, ])

# make a dataframe
gene.data2 <- data.frame(replicate = reps, time = t.n, read = y.n2)

# Plot
ggplot(gene.data2, aes(x = time, y = read)) + 
geom_point(aes(col = replicate)) +
geom_line(aes(col = replicate))



