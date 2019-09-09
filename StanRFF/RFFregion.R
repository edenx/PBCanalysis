library(randtoolbox)
library(SDALGCP)
library(sf)
library(raster)
library(tidyverse)

source("RFFfunc.R")

pop_den_ <- raster::intersect(pop_den, PBCshp) %>% replace_na(0)
PBC <- st_as_sf(PBCshp)

sf_poly_ <- raster::extract(pop_den_, PBC, cellnumbers=TRUE, small=TRUE,
                            weights=TRUE, normalizeWeights=FALSE)
lis_centr <- sapply(sf_poly_, function(x) cbind(x, scale(coordinates(pop_den_))[x[,1],]))
lis_wcentr <- sapply(lis_centr, check.0)

# check the distribution of population density 
hist(unlist(sapply(lis_wcentr, function(x) x[, "W"])))
table(unlist(sapply(lis_wcentr, function(x) sum(x[, "W"]))))
table(unlist(sapply(lis_wcentr, function(x) sum(x[, "W"]*x[,"weight"]))))

count <- PBC$X


# precompute the RFF for each lengthscale
alphas <- seq(0.2, 1, length.out = 30)
lis_Phi <- list()
lis_Phi_ <- list()
for(i in 1:length(alphas)){
        alpha <- alphas[i]
        lis_Phi[[i]] <- sim_rff(lis_wcentr, alpha=alpha)
        lis_Phi_[[i]] <- cbind(lis_Phi[[i]], PBC$Income, PBC$Crime, PBC$Environment, PBC$Employment, 
                               PBC$Barriers, PBC$propmale, PBC$Education)
}

