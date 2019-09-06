library(randtoolbox)
library(SDALGCP)
library(sf)
library(raster)
library(tidyverse)
library(tictoc)
library(rstanarm)
library(ggplot2)
library(glmnet)
library(dplyr)
library(reshape2)

source("RFFfunc.R")

        # ====================================================================== #
        #       find empirical population density distribution per polygon
        # ====================================================================== #

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


        # ====================================================================== #
        #                           find Fourier Features                        #
        # ====================================================================== #

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

        # ====================================================================== #
        #                   find lengthscle with regularised GLM                 #
        # ====================================================================== #

log_pop <- log(PBC$pop)
# the best lengthscale with regularised glm
cv_output <- alpha_cv(lis_Phi_, tune_param=0.3, offset=log_pop)
min_index <- cv_output$min_index
best_alpha <- cv_output$best_ls
best_pred <- cv_output$best_pred

hist(best_pred, main=paste0("Best Lengthscale is ", best_alpha))

        # ====================================================================== #
        #                   fit Stan model with selected lengthscale             #
        # ====================================================================== #
options(mc.cores = parallel::detectCores())
SEED <- 727282

dat_ <- as.data.frame(lis_Phi_[[min_index]]) %>%
        mutate(count=count
               , log_pop=scale(log(PBC$pop), center=FALSE)
        )

tic(paste0("Model fitting with lengthscale=", alphas[min_index]))

stan_mod_ <- stan_glm(count ~ . -log_pop,
                      offset=log_pop,
                      data=dat_, family=poisson,
                      prior=normal(0, sqrt(2)), prior_intercept=normal(0,5),
                      control = list(max_treedepth = 20),
                      # chains=5, 
                      # thin=5,
                      seed=SEED, verbose=TRUE)
toc()

yrep_ <- posterior_predict(stan_mod_)

        # ====================================================================== #
        #                   Evaluate the result with some metrics                #
        # ====================================================================== #

# PLOT:: spatial prediction with the ridge glm selected lengthscale
# plot_pred(exp(stan_mod$linear.predictors))
plot_pred(colMeans(yrep_), round(best_alpha, 2))
bias(yrep_)
rmse(yrep_)
cp95(yrep_)
