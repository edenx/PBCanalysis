library(randtoolbox)
library(SDALGCP)
library(sf)
library(raster)
library(tidyverse)
library(tictoc)
library(glmnet)

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
        phi <- sqrt(1/m) * colSums(t(w*h) %*% cbind(cos(proj/alpha), sin(proj/alpha)))

        return(phi)
}

# using the inversion -- quasi monte carlo integration;
# Matern 5-2 kernel -- student 5-2 sdf
sim_rff <- function(m=100, nu=5, alpha, plot=FALSE){
        Omega <- qt(halton(m,2), nu)
        Phi <- t(sapply(lis_wcentr, rff.region, Omega, m, alpha))
        
        Kernel <- Phi %*% t(Phi)
        if(plot){
                print(lattice::levelplot(Kernel))
        }
        
        return(Phi)
}

        # ====================================================================== #
        #                   inference with RFF: Ridge Regression
        # ====================================================================== #
# precompute the RFF for each lengthscale
alphas <- seq(0.01, 1, length.out=10)
lis_Phi <- list()
lis_Phi_ <- list()
for(i in 1:length(alphas)){
        alpha <- alphas[i]
        lis_Phi[[i]] <- sim_rff(alpha=alpha)
        lis_Phi_[[i]] <- cbind(Phi, PBC$Income, PBC$Crime, PBC$Environment, PBC$Employment, 
                      PBC$Barriers, PBC$propmale, PBC$Education)
}

# find the best lengthscale
lis_cvfit <- c()
lis_f <- list()
lis_fit <- list()

for(i in 1:length(lis_Phi)){
        Phi_ <- lis_Phi_[[i]]
        print(i)
        
        # fit a regularised glm with log link
        lambdas <- 10^seq(5, -5, length.out=100)
        cv_fit <- cv.glmnet(Phi_, count, family="poisson", 
                            # offset=log_pop,
                            alpha = 0, lambda = lambdas)
        opt_lambda <- cv_fit$lambda.min
        lis_cvfit <- c(lis_cvfit, cv_fit$cvm[which(cv_fit$lambda == opt_lambda)])
        
        cat("The optimal lambda is ", opt_lambda)
        cat("\n")
        
        fit <- glmnet(Phi_, count, family="poisson", 
                      # offset=log_pop,
                      alpha = 0, lambda = opt_lambda)
        lis_fit[[i]] <- fit
        
        lis_f[[i]] <- predict(fit, s = opt_lambda, newx = Phi_, 
                     # offset=rep(1, nrow(Phi_)),
                     type="response"
                     )
        
        hist(lis_f[[i]])
}

plot(alphas, lis_cvfit, type="line")

min_error <- which.min(lis_cvfit)
best_alpha <- alphas[min_error]
best_pred <- lis_f[[min_error]]
hist(best_pred, main=paste0("Best Lengthscale is ", best_alpha))

        # ====================================================================== #
        #                   inference with RFF: Bayesian
        # ====================================================================== #

library(rstanarm)
options(mc.cores = parallel::detectCores())

# for some length scale combine the data set
dat <- as.data.frame(lis_Phi_[8]) %>%
        mutate(count=count
               # , pop=PBC$pop
        )

SEED <- 727282

tic("Stan")
stan_glm1 <- stan_glm(count ~ ., 
                      # offset=log(pop),
                      data=dat, family=poisson, 
                      prior=normal(0, 2.5), prior_intercept=normal(0,5),
                      control = list(max_treedepth = 20),
                      # chains=5, 
                      # thin=5,
                      seed=SEED, verbose=TRUE)
toc(log=TRUE)

hist(as.vector(coef(stan_glm1)))
quantile(as.vector(coef(stan_glm1)))
hist(as.vector(se(stan_glm1)))
quantile(as.vector(se(stan_glm1)))
# predict
yrep <- posterior_predict(stan_glm1)

prop_zero <- function(y) mean(y == 0)
(prop_zero_test1 <- pp_check(stan_glm1, plotfun = "stat", stat = "prop_zero", binwidth=0.5))



