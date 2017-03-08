library(ggplot2)
library(wesanderson)
library(RColorBrewer) # display.brewer.all()
library(reshape2) # melt function
library(raster)

scale_fill_brewer(palette="Spectral")

# ggplot(faithfuld, aes(waiting, eruptions)) +
# geom_raster(aes(fill = density), interpolate = TRUE) +
# scale_fill_gradientn(colours = terrain.colors(10))
#
# ggplot(faithfuld, aes(waiting, eruptions)) +
# geom_raster(aes(fill = density), interpolate = TRUE) +
# scale_fill_gradientn(colours = rainbow(20)[1:10])
#
# pal <- wes_palette("Zissou", 11, type = "continuous")
#
# w <- ggplot(faithfuld, aes(waiting, eruptions)) +
# 	 geom_raster(aes(fill = density), interpolate = TRUE) +
# 	 scale_fill_gradientn(colours = pal)
# w

cols <- colorRampPalette(brewer.pal(10, "Spectral"))(10)
cols <- rev(cols)

v <- ggplot(faithfuld, aes(waiting, eruptions, z = density)) +
	geom_raster(aes(fill = density), interpolate = T) +
    geom_contour(colour = "white") +
	scale_fill_distiller(palette = "Spectral")
	
w <- ggplot(faithfuld, aes(waiting, eruptions, z = density)) +
	geom_raster(aes(fill = density), interpolate = T) +
    geom_contour(colour = "white") +
	scale_fill_distiller(palette = "RdYlBu")

	
	scale_fill_gradientn("Label", colours = cols)

####
# geom_point with color
####

set.seed(133)
df <- data.frame(xval=rnorm(50), yval=rnorm(50))

# Make color depend on yval
ggplot(df, aes(x=xval, y=yval, colour=yval)) + geom_point()


cols <- colorRampPalette(brewer.pal(10, "Spectral"))(20)

# Solid circles
ggplot(df, aes(x=xval, y=yval, colour = yval)) + geom_point(size = 2) + 
	scale_colour_gradientn(colours=cols) +
	theme_bw()

# Filled circles with black outline
ggplot(df, aes(x = xval, y = yval, fill = yval)) + geom_point(pch = 21, size = 5) + 
scale_fill_gradientn(colours=cols) +
theme_bw()




# Raster package
filename <- system.file("external/test.grd", package="raster")
r <- raster(filename)
r.p <- rasterToPoints(r)
df <- data.frame(r.p)

u <- ggplot(df, aes(x = x, y = y, z = test)) +
	geom_raster(aes(fill = density), interpolate = T) +
    geom_contour(colour = "white") +
	scale_fill_gradientn("Label", colours = cols)
