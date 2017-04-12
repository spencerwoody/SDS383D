###########################################################
######### Created by Spencer Woody on 27 Mar 2017 #########
###########################################################

rm(list=ls())

library(ggplot2)
library(lattice)
library(RColorBrewer) # display.brewer.all()
library(Matrix)
library(sparseMVN)
library(mvtnorm)
library(gridExtra)

col1 <- "firebrick2"
col2 <- "springgreen4"
col3 <- "dodgerblue2"

### --------------------------------------------------------------------------
### Data prep
### --------------------------------------------------------------------------

drosophila <- read.csv("Drosophila/droslong.csv", header = T)

xyplot(log2exp ~ time | gene, data = drosophila)

myfacet <- ggplot(drosophila, aes(time, log2exp)) +
geom_point(pch = 1, aes(col = replicate)) +
scale_colour_manual("Replicate", values = c(col1, col2, col3)) +
facet_wrap(~gene) 
myfacet

ggplot(drosophila, aes(time, log2exp)) +
geom_point(aes(col = replicate)) +
facet_wrap(~group)

