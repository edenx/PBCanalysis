library(randtoolbox)
library(SDALGCP)
library(sf)
library(raster)
library(tidyverse)

        # ====================================================================== #
        #       find empirical population density distribution per polygon
        # ====================================================================== #

pop_den_ <- raster::intersect(pop_den, PBCshp) %>% replace_na(0)
PBC <- st_as_sf(PBCshp)

sf_poly_ <- raster::extract(pop_den_, PBC, cellnumbers=TRUE, small=TRUE,
                            weights=TRUE, normalizeWeights=FALSE)
# may be unnecessary info
cell_index <- sapply(sf_poly_, function(x) x[,1])
poly_index <- rep(seq(1,545), times=sapply(cell_index, length))

lis_centr <- sapply(sf_poly_, function(x) cbind(x, scale(coordinates(pop_den_))[x[,1],]))

# for every intersected cells with polygon, normalise the population density per polygon
check.0 <- function(x) {
        w_val <- x[,"value"]*x[,"weight"]
        if(any(x[,"value"] != 0)){
                cbind(x, W=x[,"value"]/sum(w_val))}
        else{cbind(x, W=x[,"value"])}
}
lis_wcentr <- sapply(lis_centr, check.0)

# check the distribution of population density 
hist(unlist(sapply(lis_wcentr, function(x) x[, "W"])))
table(unlist(sapply(lis_wcentr, function(x) sum(x[, "W"]))))
table(unlist(sapply(lis_wcentr, function(x) sum(x[, "W"]*x[,"weight"]))))


# names(lis_wcentr) <- seq(1,length(lis_wcentr))

# sp_count <- raster::extract(sim_lgcp, PBC, small=TRUE, sp=TRUE, fun=sum)
# count <- sp_count[]$NZ

count <- PBC$X


        # ====================================================================== #
        #                       find Fourier Features
        # ====================================================================== #

rff.region <- function(wcentr, Omega, m, alpha){
        centr <- as.matrix(wcentr[,c("x", "y")])
        # population density empirical distribution
        w <- as.matrix(wcentr[,c("W")])
        h <- as.matrix(wcentr[,c("weight")])

        # Projection - combine data with sample frequencies
        proj <- centr %*% t(Omega) 
        
        # Fourier feature for a given area
        # print(w*h)
        phi <- sqrt(1/m) * colSums(t(w*h) %*% cbind(cos(proj/alpha), sin(proj/alpha)))

        return(phi)
}

# using the inversion -- quasi monte carlo integration;
# Matern 5-2 kernel -- student 5-2 sdf
sim_rff <- function(lis_region, m=100, nu=5, alpha, plot=FALSE){
        Omega <- qt(halton(m,2), nu)
        Phi <- t(sapply(lis_region, rff.region, Omega, m, alpha))
        
        Kernel <- Phi %*% t(Phi)
        if(plot){
                print(lattice::levelplot(Kernel))
        }
        
        return(Phi)
}

# precompute the RFF for each lengthscale
alphas <- seq(0.01, 0.3, length.out = 15)
lis_Phi <- list()
lis_Phi_ <- list()
for(i in 1:length(alphas)){
        alpha <- alphas[i]
        lis_Phi[[i]] <- sim_rff(lis_wcentr, alpha=alpha)
        lis_Phi_[[i]] <- cbind(lis_Phi[[i]], PBC$Income, PBC$Crime, PBC$Environment, PBC$Employment, 
                               PBC$Barriers, PBC$propmale, PBC$Education)
}
