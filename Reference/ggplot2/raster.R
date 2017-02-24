library(ggplot2)
library(wesanderson)
library(RColorBrewer) # display.brewer.all()
library(reshape2) # melt function
library(raster)

scale_fill_brewer(palette="Spectral")

ggplot(faithfuld, aes(waiting, eruptions)) +
geom_raster(aes(fill = density), interpolate = TRUE) +
scale_fill_gradientn(colours = terrain.colors(10))

ggplot(faithfuld, aes(waiting, eruptions)) +
geom_raster(aes(fill = density), interpolate = TRUE) +
scale_fill_gradientn(colours = rainbow(20)[1:10])

pal <- wes_palette("Zissou", 11, type = "continuous")

w <- ggplot(faithfuld, aes(waiting, eruptions)) +
	 geom_raster(aes(fill = density), interpolate = TRUE) +
	 scale_fill_gradientn(colours = pal)
w

cols <- colorRampPalette(brewer.pal(8,"Dark2"))(length(table.data))

cols <- colorRampPalette(brewer.pal(10, "Spectral"))(10)
cols <- rev(cols)


v <- ggplot(faithfuld, aes(waiting, eruptions, z = density)) +
	geom_raster(aes(fill = density), interpolate = T) +
    geom_contour(colour = "white") +
	scale_fill_distiller(palette = "Spectral")
	
	
	scale_fill_gradientn("Label", colours = cols)



something.m <- melt(something)

p <- ggplot(something.m, aes(x = Var2, y = Var1, fill = value)) + # change fill to z?
geom_tile() + 
scale_fill_gradient(low = "snow", high = "dodgerblue3") + 
xlab("Time point") + 
ylab("Cluster number") +
ggtitle("RNA Cluster Centroids")
p

scale_colour_brewer(palette = "Greens")


# Raster package
filename <- system.file("external/test.grd", package="raster")
r <- raster(filename)
r.p <- rasterToPoints(r)
df <- data.frame(r.p)

u <- ggplot(df, aes(x = x, y = y, z = test)) +
	geom_raster(aes(fill = density), interpolate = T) +
    geom_contour(colour = "white") +
	scale_fill_gradientn("Label", colours = cols)
