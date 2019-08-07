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


names(lis_wcentr) <- seq(1,length(lis_wcentr))

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
        phi <- sqrt(1/m) * colSums(t(w*h) %*% cbind(cos(proj), sin(proj)))
        # phi <- sqrt(1/m) * cbind(cos(proj/alpha), sin(proj/alpha))

        return(phi)
}

# using the inversion -- quasi monte carlo integration;
# Matern 5-2 kernel -- student 5-2 sdf
alpha <- 0.25   # lengthscale
m <- 1000
Omega <- qt(halton(m,2), 5)/alpha

Phi <- t(sapply(lis_wcentr, rff.region, Omega, m, alpha))
hist(Phi)
corr <- cor(Phi)
dev.new(width = 1000, height = 700, unit = "px")
iplotCorr(corr)
lattice::levelplot(corr)

Kernel <- Phi %*% t(Phi) # approximation of Kernel matrix (regional level)
Kernel[1:5, 1:5]
lattice::levelplot(Kernel)
        # ====================================================================== #
        #                   inference with RFF: Ridge Regression
        # ====================================================================== #

Phi_ <- cbind(Phi, PBC$Income, PBC$Crime, PBC$Environment)
# fit a regularised glm with log link
library(glmnet)
lambdas <- 10^seq(5, -5, length.out=100)
cv_fit <- cv.glmnet(Phi, count, family="poisson", alpha = 0, lambda = lambdas)
opt_lambda <- cv_fit$lambda.min
opt_lambda

fit <- glmnet(Phi_, count, family="poisson", alpha = 0, lambda = opt_lambda)
f <- predict(fit, s = opt_lambda, newx = Phi_, type="response")

# clarely there's overdispersion; 
# should account for this using Negative Binomial with log link

        # ====================================================================== #
        #                   inference with RFF: Bayesian
        # ====================================================================== #



