library(sf)
library(SDALGCP)
library(raster)
library(tibble)
library(dplyr)
library(tidyverse)
library(ggplot2)

# extracting raster object for the Newcastle upon tyne region of interest
pop_den_ <- raster::intersect(pop_den, PBCshp) %>% replace_na(0)
plot(pop_den_, axes=TRUE)
head(pop_den_[])

# transform `sp` to `sf`
PBC <- st_as_sf(PBCshp)
# notably, population density is a measurement of population per unit area

pop_ <- ggplot() + 
        geom_sf(data=PBC["pop"], aes(fill=pop), lwd=0.05) + 
        scale_fill_distiller(palette="BuPu", direction=1) +
        theme_linedraw() +
        labs(x="Northing", y="Easting")

count_ <- ggplot() + 
        geom_sf(data=PBC["X"], aes(fill=X), lwd=0.05) + 
        scale_fill_distiller(palette="BuPu", direction=1) +
        theme_linedraw() +
        labs(x="Northing", y="Easting")

# create data frame for computational grid
plot_rast <- as_tibble(coordinates(pop_den_)) %>% mutate(den=pop_den_[]) 
rast_ <- ggplot(plot_rast) + geom_raster(aes(x=x, y=y, fill=den), interpolate=TRUE) + 
        scale_fill_distiller(palette="BuPu", direction=1) +
        geom_sf(data=PBC["geometry"], fill=NA, lwd=0.1) +
        theme_linedraw() +
        labs(x="Northing", y="Easting")
# combine plots
gridExtra::grid.arrange(pop_, count_, rast_, nrow=3)


