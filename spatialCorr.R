library(rspatial)
library(geodata)
library(sp)
library(spData)

datapath <- "X:/georisk/HaRIA_B_Wind/data/derived/obs/1-minute/wind"
plotpath <- "X:/georisk/HaRIA_B_Wind/data/derived/obs/1-minute/wind/gammagpdcon"
paramdf <- read.csv(paste(plotpath, "params.csv", sep = "/"))

dsp <- SpatialPoints(paramdf[,3:2], proj4string = CRS("+proj=longlat +datum=GDA2020"))

dsp <- SpatialPointsDataFrame(dsp, paramdf)
blues <- colorRampPalette(c('yellow', 'orange', 'green', 'blue', 'dark blue'))
cuts <- c(0,5,10,15,20,25,30)
xicuts <- c(-1., -0.5, -0.25, -0.1, -0.05, -0.01, 0, 0.5, 1.0, 2.)
aus <- gadm("AUS", level=1, path=tempdir())
pols <- list("sp.polygons", aus, fill='lightgray')
spplot(dsp, 'xi', cuts=xicuts,  sp.layout=pols, pch=20, cex=2)
