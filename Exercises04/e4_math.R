###########################################################
######### Created by Spencer Woody on 11 Feb 2017 #########
###########################################################

library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer) # display.brewer.all()
library(xtable)

# Read in the data

math <- read.csv("Math/mathtest.csv", header = T)

# Boxplot

box <- ggplot(math, aes(factor(school), mathscore)) +
geom_boxplot()

# Sample size vs. group average

math.df <- data.frame(cbind(
school = unique(math$school),
count = as.vector(table(math$school)),
score.av = as.vector(tapply(math$mathscore, math$school, mean)),
score.sd = as.vector(tapply(math$mathscore, math$school, sd))
))

head(math.df)

scatter <- ggplot(math.df, aes(count, score.av)) + 
geom_point() +
ggtitle("Average Test Scores vs. Sample Size for each school") +
xlab("Sample size") +
ylab("Average Math Test Score") +
theme_bw() +
theme(plot.title = element_text(hjust = 0.5)) 

pdf("Math/img/scatter.pdf")
scatter
dev.off()

## ---------------------------------------------------------------------------
## Gibbs sampler 
## ---------------------------------------------------------------------------

I <- nrow(math.df)
y <- math$mathscore
Ni <- math.df$count
N <- sum(Ni)
ybar <- math.df$score.av
Ni.ybar <- Ni * ybar

niter <- 1e4
nburn <- 2e3

thetas <- matrix(nrow = niter, ncol = I)
lambda <- rep(NA, niter)
gamma  <- rep(NA, niter)
mu     <- rep(NA, niter)

# Initialize Gibbs sampler
thetas[1, ] <- rep(0, I)
lambda[1]   <- 1
gamma[1]    <- 1
mu[1]       <- 0

for (iter in 2:niter) {
	# Update thetas
	post.var  <- (Ni * lambda[iter - 1] + 
			     lambda[iter - 1] * gamma[iter - 1]) ^ -1
	post.mean <- post.var * (lambda[iter - 1] * Ni.ybar + 
							 lambda[iter - 1] * gamma[iter - 1] * mu[iter - 1])
	
	thetas[iter, ] <- rnorm(I, post.mean, sqrt(post.var))
	
	# Update mu
	mu[iter] <- rnorm(1, 
				mean(thetas[iter, ]),
				(lambda[iter - 1] * gamma[iter - 1] * I) ^ -0.5)
					
	# Update lambda
	S.y <- sum( (y - rep(thetas[iter, ], times = Ni)) ^ 2 )
	S.theta <- sum( (thetas[iter, ] - mu[iter]) ^ 2 )
	rate.lambda <- (S.y + gamma[iter - 1] * S.theta) / 2
	
	lambda[iter] <- rgamma(1, 
				     	   (N + I) / 2,
						   rate.lambda)
						   
	# Update gamma
	rate.gamma <- lambda[iter] * sum((thetas[iter, ] - mu[iter]) ^ 2) / 2
	
	gamma[iter] <- rgamma(1, 
						  I / 2 - 1,
						  rate.gamma)
						  
	if (iter %% 1000 == 0) {
		print(sprintf("Iteration %i out of %i...", iter, niter))
	}
}

# Remove burn-in
thetas <- thetas[(nburn+1):niter, ]
mu <- mu[(nburn+1):niter]
lambda <- lambda[(nburn+1):niter]
gamma <- gamma[(nburn+1):niter]

# plot(mu, type = "l")
# plot(lambda, type = "l")
# plot(gamma, type = "l")

# 95% credible intervals
post.summary <- rbind(
	quantile(mu, c(0.025, 0.5, 0.975)),
	quantile(lambda, c(0.025, 0.5, 0.975)),
	quantile(gamma, c(0.025, 0.5, 0.975))
	)
rownames(post.summary) <- c("mu", "lambda", "gamma")
post.summary

# Make plots of posterior distributions
mu.plot <- qplot(mu, bins = 15, fill = I("grey30"), col = I("grey90")) + ggtitle("Posterior distribution of mu") + 
xlab("mu") +
theme_bw() +
theme(plot.title = element_text(hjust = 0.5)) 

lambda.plot <- qplot(lambda, bins = 15, fill = I("grey30"), col = I("grey90")) + ggtitle("Posterior distribution of lambda") + 
xlab("lambda") +
theme_bw() +
theme(plot.title = element_text(hjust = 0.5)) 

gamma.plot <- qplot(gamma, bins = 15, fill = I("grey30"), col = I("grey90")) + ggtitle("Posterior distribution of gamma") + 
xlab("gamma") +
theme_bw() +
theme(plot.title = element_text(hjust = 0.5)) 

# Plot shrinkage coefficients
thetas.hat <- apply(thetas, 2, mean)
kappa <- abs(1 - thetas.hat / ybar)

math.df2 <- data.frame(cbind(math.df, kappa = kappa))

kappa.plot <- ggplot(math.df2, aes(count, kappa)) + 
geom_point() +
ggtitle("Shrinkage coefficient vs. Sample Size for each school") +
xlab("Sample size") +
ylab("Absolute shrinkage coefficient") +
theme_bw() +
theme(plot.title = element_text(hjust = 0.5)) 

pdf("Math/img/kappa.pdf")
kappa.plot
dev.off()

