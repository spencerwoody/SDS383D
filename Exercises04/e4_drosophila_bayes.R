###########################################################
######### Created by Spencer Woody on 27 Mar 2017 #########
###########################################################

# To-do list:
#    1. Tidy up the code (document functions and )
#    2. Function for fitting and plotting
#    3. Function for looping over all clusters, genes, and replicates


rm(list=ls())

library(ggplot2)
library(lattice)
library(RColorBrewer) # display.brewer.all()
library(Matrix)
library(mvtnorm)
library(gridExtra)

source("e4_covmat_funs.R")

### --------------------------------------------------------------------------
### Data prep
### --------------------------------------------------------------------------

drosophila <- read.csv("Drosophila/droslong.csv", header = T)

# Simpler dataframe
fruitfly <- data.frame(
	group = drosophila$group, 
	gene = drosophila$gene,
	replicate = drosophila$replicate, 
	time = drosophila$time,
	log2exp = drosophila$log2exp)
	
# Order the dataframe
fruitfly <- fruitfly[with(fruitfly, order(group, gene, replicate, time)), ]

head(fruitfly)

### --------------------------------------------------------------------------
### Optimize parameters
### --------------------------------------------------------------------------

# Choose a group / cluster
groupnum <- "1"
mygroup <- paste0("group", groupnum)

# Choose all data in this cluster
Yi <- fruitfly[fruitfly$group == mygroup, ] #
genes.i <- unique(Yi$gene)
Ni <- length(unique(Yi$gene)) #

# Response vector
yvec <- Yi$log2exp

# Initial parameters (on log scale)
gamma.f <- -2 # IN exp
alpha.f <- -2 # OUT exp
gamma.g <- 1
alpha.g <- 1
gamma.h <- 2
alpha.h <- 2
sigma2 <- -2

init.par <- c(gamma.f, alpha.f, gamma.g, alpha.g, gamma.h, alpha.h, sigma2)

# Choose one gene in cluster dataframe to make
# vector of times for samples, and number of samples for each replicate
yn <- Yi[Yi$gene == Yi$gene[1], ] #

reps <- unique(Yi$replicate) #
numreps <- length(reps) #

tnr.list <- list()
Nnr <- NA

for (k in 1:numreps) {
	tnr.list[[ k ]] <- yn$time[yn$replicate == reps[k]]
	Nnr[k] <- length(tnr.list[[ k ]])
}

tn <- unlist(tnr.list) # Concatenated times over all replicates
D <- length(tn) # Number of measurements across all replicates

# Find optimal parameters by maximizing the log-likelihood
myopt <- optim(par = init.par, fn = neg.ll, Yi = Yi)
myopt

optpar <- myopt$par 

# For group 1... the log-likelihood as a function of the 
# covariance function parameters is very flat, so alpha.h can shoot up
if (optpar[6] > 10) {
	optpar[6] <- 10
}

optCOV <- full.covmat(Yi, optpar)
optCOV.inv <- solve(optCOV)

# Sequence of times to smooth over
t.star <- seq(0.5, 12.5, length.out = 200)

### --------------------------------------------------------------------------
### Cluster-level function (denoted by i)
### --------------------------------------------------------------------------

# Parameters for cluster-level function
params.i <- c(exp(optpar[c(5,6)]), 0)

# Off-diagonal block of covariance matrix of data and g.star
Ki.star <- make.covmat(t.star, rep(tn, Ni), cov.fun = C.SE, params = params.i)
# Marginal covariance of g.star
Ki.starstar <- make.covmat(t.star, cov.fun = C.SE, params = params.i)

# Posterior mean of h function
hi.est <- Ki.star %*% optCOV.inv %*% yvec

# Posterior covariance matrix of hi.est
hi.covmat <- Ki.starstar - Ki.star %*% optCOV.inv %*% t(Ki.star)
hi.var <- diag(hi.covmat)

# 95% confidence intervals
hi.lo <- hi.est - 1.96 * sqrt(hi.var)
hi.hi <- hi.est + 1.96 * sqrt(hi.var)

# Create a dataframe for plotting
hi.df <- data.frame(time = t.star, lo = hi.lo, md = hi.est, hi = hi.hi)

# Plot
hi.plot <- ggplot(hi.df, aes(time)) + 
geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70", alpha = 0.5) +
geom_line(aes(y = md)) +
geom_point(data = drosophila[drosophila$group == mygroup, ], aes(time, log2exp, col = gene)) +
scale_colour_brewer(palette="Set1") +
scale_x_continuous(breaks = seq(0, 12, by = 2)) +
ggtitle(paste0("Estimation of Group level Function for Group ", groupnum)) +
xlab("Time") +
ylab("log2-Normalized gene expression") + 
theme_bw() +
theme(plot.title = element_text(hjust = 0.5)) 
hi.plot

### --------------------------------------------------------------------------
### Gene- and replicate-level function
### --------------------------------------------------------------------------

gn.df.list <- list()
fr.df.list <- list()

for (k in 1:Ni) {
	# Which part of off-diagonal covariance matrix must be changed?
	gene_indices <- ((k - 1) * D + 1):(k * D)

	# Gene-level covariance function parameters
	params.n <- c(exp(optpar[3:4]), 0)

	# Make gene-level off-diagonal covariance matrix
	Kn.starsmall <- make.covmat(t.star, tn, cov.fun = C.SE, params = params.n)
	Kn.star <- Ki.star
	Kn.star[, gene_indices] <- Kn.star[, gene_indices] + Kn.starsmall

	# Marginal covariance of gn.star, 
	Kn.starstar <- Ki.starstar + 
				   make.covmat(t.star, cov.fun = C.SE, params = params.n)

	# Posterior mean of gn function
	gn.est <- Kn.star %*% optCOV.inv %*% yvec

	# Conditional covariance of gn
	gn.cov <- Kn.starstar - Kn.star %*% optCOV.inv %*% t(Kn.star)
	gn.var <- diag(gn.cov)

	# 95% confidence intervals
	gn.lo <- gn.est - 1.96 * sqrt(gn.var)
	gn.hi <- gn.est + 1.96 * sqrt(gn.var)

	gn.df.k <- data.frame(
		gene = genes.i[k], 
		time = t.star, 
		lo = gn.lo, 
		est = gn.est, 
		hi = gn.hi)
	
	gn.df.list[[ k ]] <- gn.df.k
	
	for (l in 1:numreps) {
		tnr <- tnr.list[[ l ]]
		
		params.r <- c(exp(optpar[1:2]), 0)
		
		if (l == 1) {
			rep_indices <- 1:Nnr[l]
		} else if (l == numreps) {
			rep_indices <- tail(1:D, n = Nnr[l])
		} else {
			sum_l <- sum(Nnr[which(1:numreps < l)])
			rep_indices <- (sum_l + 1):(sum_l + Nnr[l])
		}
		
		
		# Make gene-level off-diagonal covariance matrix
		Kr.starsmall <- make.covmat(t.star, tnr, 
								    cov.fun = C.SE, params = params.r)
		Kr.star <- Kn.star
		Kr.star[, gene_indices[rep_indices]] <- 
			Kr.star[, gene_indices[rep_indices]] + Kr.starsmall
		
		Kr.starstar <- Kn.starstar + 
					   make.covmat(t.star, cov.fun = C.SE, params = params.r)
			
		# Posterior mean of fnr function		   
		fr.est <- Kr.star %*% optCOV.inv %*% yvec
		
		# Conditional covariance of fnr
		fr.cov <- Kr.starstar - Kr.star %*% optCOV.inv %*% t(Kr.star)
		fr.var <- diag(fr.cov)

		# 95% confidence intervals
		fr.lo <- fr.est - 1.96 * sqrt(fr.var)
		fr.hi <- fr.est + 1.96 * sqrt(fr.var)

		fr.df.kl <- data.frame(
			replicate = reps[l],
			gene = genes.i[k], 
			time = t.star, 
			lo = fr.lo, 
			est = fr.est, 
			hi = fr.hi)
	
		fr.df.list[[ length(fr.df.list) + 1 ]] <- fr.df.kl
	}
}

gn.df.full <- do.call("rbind", gn.df.list)
fr.df.full <- do.call("rbind", fr.df.list)

gn.plot.full <- ggplot(gn.df.full, aes(time)) + 
geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70", alpha = 0.8) +
geom_line(aes(y = est)) +
geom_point(data = drosophila[drosophila$group == mygroup, ], aes(time, log2exp, col = replicate)) +
facet_wrap(~gene) + 
scale_colour_brewer(palette="Set1") +
scale_x_continuous(breaks = seq(0, 12, by = 2)) +
ggtitle(paste0("Estimation of Gene level Functions for Group ", groupnum)) +
xlab("Time") +
ylab("log2-Normalized gene expression") + 
theme_bw() +
theme(plot.title = element_text(hjust = 0.5)) 
gn.plot.full

fr.plot.full <- ggplot(fr.df.full, aes(time)) + 
geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70", alpha = 0.8) +
geom_line(aes(y = est)) +
geom_point(data = drosophila[drosophila$group == mygroup, ], aes(time, log2exp, col = gene)) +
facet_grid(gene ~ replicate) + 
scale_colour_brewer(palette="Set1") +
scale_x_continuous(breaks = seq(0, 12, by = 2)) +
ggtitle(paste0("Estimation of Gene-Replicate level Functions for Group ", groupnum)) +
xlab("Time") +
ylab("log2-Normalized gene expression") + 
theme_bw() +
theme(plot.title = element_text(hjust = 0.5)) 
fr.plot.full

pdf(paste0("Drosophila/img/GPgroup", groupnum, "_group.pdf"), 
width = 7.5, height = 7)
hi.plot
dev.off()

pdf(paste0("Drosophila/img/GPgroup", groupnum, "_genes.pdf"), 
width = 9, height = 7)
gn.plot.full
dev.off()

pdf(paste0("Drosophila/img/GPgroup", groupnum, "_replicates.pdf"), 
width = 6, height = 1.6 * Ni)
fr.plot.full
dev.off()


### --------------------------------------------------------------------------
### Plot 
### --------------------------------------------------------------------------



# plot

allgenes <- ggplot(drosophila, aes(time, log2exp)) +
geom_point(aes(col = replicate)) +
# scale_colour_manual("Replicate", values = c(col1, col2, col3)) +
facet_wrap(~gene) +
ggtitle("Expression Profiles for all Genes") +
xlab("Time") +
ylab("log2-Normalized gene expression") + 
scale_colour_brewer(palette="Set1") +
scale_x_continuous(breaks = seq(0, 12, by = 2)) +
theme_bw() +
theme(plot.title = element_text(hjust = 0.5))
allgenes

pdf("Drosophila/img/allgenes.pdf")
allgenes
dev.off()

# Facet plots for each group

myfacet1 <- ggplot(drosophila[drosophila$group == "group1", ], aes(time, log2exp)) +
geom_point(aes(col = replicate)) +
# scale_colour_manual("Replicate", values = c(col1, col2, col3)) +
facet_wrap(~gene) +
ggtitle("Expression Profiles for Genes in Group 1") +
xlab("Time") +
ylab("log2-Normalized gene expression") + 
scale_colour_brewer(palette="Set1") +
scale_x_continuous(breaks = seq(0, 12, by = 2)) +
theme_bw() +
theme(plot.title = element_text(hjust = 0.5))
myfacet1

pdf("Drosophila/img/facetgroup1.pdf", width = 5.5, height = 5)
myfacet1
dev.off()

myfacet2 <- ggplot(drosophila[drosophila$group == mygroup, ], aes(time, log2exp)) +
geom_point(aes(col = replicate)) +
facet_wrap(~gene) +
ggtitle("Expression Profiles for Genes in Group 2") +
xlab("Time") +
ylab("log2-Normalized gene expression") + 
scale_colour_brewer(palette="Set1") +
scale_x_continuous(breaks = seq(0, 12, by = 2)) +
theme_bw() +
theme(plot.title = element_text(hjust = 0.5))
myfacet2

pdf("Drosophila/img/facetgroup2.pdf", width = 5.5, height = 5)
myfacet2
dev.off()


myfacet3 <- ggplot(drosophila[drosophila$group == "group3", ], aes(time, log2exp)) +
geom_point(aes(col = replicate)) +
facet_wrap(~gene) +
ggtitle("Expression Profiles for Genes in Group 3") +
xlab("Time") +
ylab("log2-Normalized gene expression") + 
scale_colour_brewer(palette="Set1") +
scale_x_continuous(breaks = seq(0, 12, by = 2)) +
theme_bw() +
theme(plot.title = element_text(hjust = 0.5))
myfacet3

pdf("Drosophila/img/facetgroup3.pdf", width = 5.5, height = 5)
myfacet3
dev.off()

# All genes, faceted by cluster

allgenefacet <- ggplot(drosophila, aes(time, log2exp)) +
geom_point(aes(col = replicate)) +
ggtitle("Expression Profiles for all Genes in Group1") +
xlab("Time") +
ylab("log2-Normalized gene expression") + 
scale_colour_brewer(palette="Set1") +
scale_x_continuous(breaks = seq(0, 12, by = 2)) +
theme_bw() +
theme(plot.title = element_text(hjust = 0.5)) + 
facet_wrap(~group)

pdf("Drosophila/img/allgenefacet.pdf", width = 8, height = 3)
allgenefacet
dev.off()


allgenes1 <- ggplot(drosophila[drosophila$group == "group1", ], aes(time, log2exp)) +
geom_point(aes(col = gene)) +
ggtitle("Expression Profiles for all Genes in Group 1") +
xlab("Time") +
ylab("log2-Normalized gene expression") + 
scale_colour_brewer(palette="Set1") +
scale_x_continuous(breaks = seq(0, 12, by = 2)) +
theme_bw() +
theme(plot.title = element_text(hjust = 0.5))

pdf("Drosophila/img/allgenes1.pdf", width = 5.5, height = 5)
allgenes1
dev.off()

allgenes2 <- ggplot(drosophila[drosophila$group == mygroup, ], aes(time, log2exp)) +
geom_point(aes(col = gene)) +
ggtitle("Expression Profiles for all Genes in Group 2") +
xlab("Time") +
ylab("log2-Normalized gene expression") + 
scale_colour_brewer(palette="Set1") +
scale_x_continuous(breaks = seq(0, 12, by = 2)) +
theme_bw() +
theme(plot.title = element_text(hjust = 0.5))

pdf("Drosophila/img/allgenes2.pdf", width = 5.5, height = 5)
allgenes2
dev.off()

allgenes3 <- ggplot(drosophila[drosophila$group == "group3", ], aes(time, log2exp)) +
geom_point(aes(col = gene)) +
ggtitle("Expression Profiles for all Genes in Group 3") +
xlab("Time") +
ylab("log2-Normalized gene expression") + 
scale_colour_brewer(palette="Set1") +
scale_x_continuous(breaks = seq(0, 12, by = 2)) +
theme_bw() +
theme(plot.title = element_text(hjust = 0.5))

pdf("Drosophila/img/allgenes3.pdf", width = 5.5, height = 5)
allgenes3
dev.off()

pdf("Drosophila/img/allgenefacet.pdf", width = 8, height = 3)
allgenefacet
dev.off()


