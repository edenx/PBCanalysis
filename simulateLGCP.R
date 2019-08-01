# load required packages
library(sf)
library(SDALGCP)
library(raster)
library(tibble)
library(dplyr)
library(tidyverse)
library(ggplot2)

library(spatstat)
library(maptools)
library(RandomFields)

library(rgeos)
library(rgdal)
library(fasterize)

source("polyOpt.R")

# find polygon and pixel level population density statistics
pop_den_ <- raster::intersect(pop_den, PBCshp) %>% replace_na(0)
PBC <- st_as_sf(PBCshp)
sf_pop_den_ <- pop_den_

sf_out <- sptpolyOpt(sf_pop_den_, PBC, plot=FALSE)
sf_pixel <- sf_out[[1]]
sf_pop_den_ <- as(as(sf_pop_den_, "SpatialGridDataFrame"), "sf")

index <- as.integer(rownames(sf_pop_den_[sf_pixel$geometry,"NZ"]))
popwNZ <- rep(0, nrow(sf_pop_den_))
popwNZ[index] <- sf_pixel$normNZ
        
sf_pop_den_ <- sf_pop_den_ %>% mutate(popwNZ=popwNZ)

win <- as.owin(as(pop_den_, "SpatialGridDataFrame"))

# Data Simulation on the computational grid
model <- RMexp(var=1, scale=0.2)

w <- as.mask(w=win)
xcol <- w$xcol
yrow <- w$yrow
dimw <- w$dim
Lambda <- as.im(w)

set.seed(241)
spc <- RandomFields::RFoptions()$general$spConform
if(spc) RandomFields::RFoptions(spConform=FALSE)
sim <- RandomFields::RFsimulate(model, xcol, yrow, grid=TRUE, n=1)
if(spc) RandomFields::RFoptions(spConform=TRUE)

Lambda[] <- exp(0.1*rnorm(nrow(sf_pop_den_)) + sim[,]) * as.vector(sf_pop_den_$popwNZ)

# Given intensity, sample from homogenous Poisson for each cell
X <- sapply(Lambda$v, function(x) rpois(1, x))

# Create a raster layer to store the simulation
sim_lgcp <- pop_den_
sim_lgcp[] <- X

# visualisation
plot_sim <- as_tibble(coordinates(sim_lgcp)) %>% mutate(den=sim_lgcp[]) 
sim_ <- ggplot(plot_sim) + geom_raster(aes(x=x, y=y, fill=den), interpolate=TRUE) + 
        scale_fill_distiller(palette="BuPu", direction=1) +
        geom_sf(data=PBC["geometry"], fill=NA, lwd=0.1) +
        theme_linedraw() +
        labs(x="Northing", y="Easting")
sim_
# =======================================================================================
# aggregate to regions
# getting values from raster within polygons
# set `normalizeWeights=TRUE` to get for each cell the proportion within a polygon
extr_sim_lgcp <- raster::extract(sim_lgcp, PBC, weights=TRUE, normalizeWeights=TRUE)

smth_sim_lgcp <- sapply(extr_sim_lgcp, weighted.sim_lgcp)
sum_sim_lgcp <- sapply(extr_sim_lgcp, sum.sim_lgcp)

# assign to a new spatial polygon data frame
simPBC <- PBC[, c("X", "geometry")] %>% mutate(X=smth_sim_lgcp)


