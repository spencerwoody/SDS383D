###########################################################
######### Created by Spencer Woody on 27 Mar 2017 #########
###########################################################

library(ggplot2)
library(RColorBrewer) # display.brewer.all()
library(lme4)
library(Matrix)
library(sparseMVN)
library(mvtnorm)
library(gridExtra)

### --------------------------------------------------------------------------
### Data prep
### --------------------------------------------------------------------------

polls <- read.csv("Polls/polls.csv", header = T)

states <- polls$state
  
N <- nrow(polls)

# Design matrix
X <- cbind(rep(1, N),
polls$edu == "HS", polls$edu == "SomeColl", polls$edu == "Bacc",
polls$age == "30to44", polls$age == "45to64", polls$age == "65plus",
polls$female,
polls$black)

# Response vector
y <- polls$bush

# Create new

### --------------------------------------------------------------------------
### Hierarchical model with lme4 package
### --------------------------------------------------------------------------

model1 <- glmer(bush ~ edu + age + female + black + 
   (1 | state), data = polls, family = binomial)
 
model2 <- glm(bush ~ edu + age + female + black, data = polls, family = binomial)

