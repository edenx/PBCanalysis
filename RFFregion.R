# Given Region i, we find centroid of the computational cell x (vector length n)
n <- sample(seq(from=5, to=30), size=504, replace=TRUE)      # number of cells within each region (504 in total)
w <- sapply(n, function(x) runif(x))
lis_centr <- sapply(n, function(x) cbind(scale(seq(-2,2,length.out=x)), scale(seq(-2,2,length.out=x))))
lis_wcentr <- sapply(seq(504), function(x) cbind(lis_centr[[x]], w[[x]]))


library(randtoolbox)
library(SDALGCP)
library(sf)
library(raster)
library(tidyverse)

source("polyOpt.R")

        # ====================================================================== #
        #                        find coordinates and weights
        # ====================================================================== #

pop_den_ <- raster::intersect(pop_den, PBCshp) %>% replace_na(0)
PBC <- st_as_sf(PBCshp)
sf_pop_den_ <- pop_den_

sf_out <- sptpolyOpt(sf_pop_den_, PBC, plot=FALSE)
sf_pixel <- sf_out$sf.pixel
sf_centr <- sf_pixel[c("LSOA04CD", "X", "normNZ")]

df_centr <- as.data.frame(sf_centr) %>% 
        rename(count=X, W=normNZ) %>%
        # centering the coordinates
        cbind(., as.data.frame(scale(st_coordinates(sf_centr)))) %>%
        mutate(geometry=NULL)

# drop unused levels
                                      # WTF????? NA????
                                      # unique(sf_out$sf.poly$LSOA04CD)
                                      # E01008430 is NA!!!!!!!
                                      # But unique levels of df_centr is 544 instead of 545

# temp <- levels(droplevels(df_centr$LSOA04CD))
df_centr$LSOA04CD <- droplevels(df_centr$LSOA04CD)

# create list of pixels within each polygon (504 in total)
lis_wcentr <- split(df_centr, df_centr$LSOA04CD)
count <- sapply(split(df_centr[,2], df_centr$LSOA04CD), unique)

        # ====================================================================== #
        #                       find Fourier Features
        # ====================================================================== #

# using the inversion -- quasi monte carlo integration;
# Matern 5-2 kernel -- student 5-2 sdf
alpha <- 0.5    # bandwidth
m <- 100
                                        # Parameterisation of Matern 5-2???????????? 
Omega <- qt(halton(m,2), 5, 1/alpha)

rff.region <- function(wcentr, Omega, m, alpha){
        centr <- as.matrix(wcentr[,c("X", "Y")])
        w <- as.matrix(wcentr[,c("W")])
        
        # Projection - combine data with sample frequencies
        proj <- centr %*% t(Omega) 
        # Fourier feature for a given area
        phi <- sqrt(1/m) * colSums(t(w) %*% cbind(cos(proj/alpha), sin(proj/alpha))) 
        
        return(phi)
}

Phi <- t(sapply(lis_wcentr, rff.region, Omega, m, alpha))
Kernel <- Phi %*% t(Phi) # approximation of Kernel matrix (regional level)

        # ====================================================================== #
        #                   inference with RFF (regularised)
        # ====================================================================== #

Phi_ <- cbind(Phi, PBC$Income, PBC$Crime, PBC$Environment)
# fit a regularised glm with log link
library(glmnet)
lambdas <- 10^seq(5, -5, length.out=100)
cv_fit <- cv.glmnet(Phi, count, family="poisson", alpha = 0, lambda = lambdas)
opt_lambda <- cv_fit$lambda.min

fit <- glmnet(Phi, count, family="poisson", alpha = 0, lambda = opt_lambda)
f <- predict(fit, s = opt_lambda, newx = Phi)
