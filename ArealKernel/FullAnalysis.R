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

source("RegKernelFunc.r")

        # ====================================================================== #
        #       find empirical population density distribution per polygon
        # ====================================================================== #

pop_den_ <- raster::intersect(pop_den, PBCshp) %>% replace_na(0)
PBC <- st_as_sf(PBCshp)
PBC_df <- PBC %>% st_set_geometry(NULL)

# value: population density per cell
# weight: the proportion of cell intersected with the polygon
# W: weighted population density per grid intersection with polygon
sf_poly_ <- raster::extract(pop_den_, PBC, cellnumbers=TRUE, small=TRUE,
                            weights=TRUE, normalizeWeights=FALSE)
lis_centr <- sapply(sf_poly_, function(x) cbind(x, scale(coordinates(pop_den_))[x[,1],]))
lis_wcentr <- sapply(lis_centr, check.0)

# find the exposure risk
pop_den_reg <- sapply(lis_centr, function(x) sum(x[,"value"] * x[,"weight"]))
exposure <- PBC_df$pop/pop_den_reg

# check the distribution of population density 
hist(unlist(sapply(lis_wcentr, function(x) x[, "W"])))
table(unlist(sapply(lis_wcentr, function(x) sum(x[, "W"]))))
table(unlist(sapply(lis_wcentr, function(x) sum(x[, "W"]*x[,"weight"]))))

count <- PBC$X

        # ====================================================================== #
        #               if want to look at the exact KME                         #
        # ====================================================================== #
# choose some lengthscale alpha=0.4
# paralleled but still quite slow to compute (as a double for loop is required)
# upshot: just use RFF with 100 features
tic()
regker <- sim.regker(lis_wcentr, alpha=0.4, plot=TRUE)
toc()

tic()
inv_regker <- Matrix::solve(regker)
toc()

lattice::levelplot(inv_regker)

tic()
rffker <- sim.rff(lis_wcentr, alpha=0.4, plot=TRUE)
toc()

        # ====================================================================== #
        #                           find Fourier Features                        #
        # ====================================================================== #

# precompute the RFF for each lengthscale
alphas <- seq(0.2, 1, length.out = 30)
lis_Phi.null <- list()
lis_Phi <- list()
lis_Phi_ <- list()

for(i in 1:length(alphas)){
        alpha <- alphas[i]
        
        lis_Phi.null[[i]] <- sim.rff(lis_wcentr, alpha=alpha)
        
        lis_Phi[[i]] <- lis_Phi.null[[i]] %>% as.data.frame() %>%
                dplyr::mutate(count=PBC_df$X)
        
        lis_Phi_[[i]] <- lis_Phi[[i]] %>% 
                cbind(dplyr::select(PBC_df, Income, Crime, Environment, Employment, 
                                    Barriers, propmale, Education))
}

        # ====================================================================== #
        #                   find lengthscle with regularised GLM                 #
        # ====================================================================== #

# the best lengthscale with regularised glm
cv_output <- alpha_ridge(lis_Phi_, exposure)
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
        mutate(count=count)

tic(paste0("Model fitting with lengthscale=", alphas[min_index]))


stan_mod_ <- stan_glm(count ~ . ,
                      offset=exposure,
                      data=dat_, family=poisson,
                      prior=normal(0, sqrt(2)), prior_intercept=normal(0,5),
                      control = list(max_treedepth = 20),
                      # chains=5, 
                      # thin=5,
                      seed=SEED, verbose=TRUE)
toc()

yrep_ <- posterior_predict(stan_mod_)


        # ====================================================================== #
        #             fit INLA for the list of lengthscales (default prior)      #
        # ====================================================================== #

lmod <- lm(count ~ . , lis_Phi_[[1]])
form <- formula(lmod)
## extreme values for small ls 
## Perhaps prior with smaller variance?
## (so that extrme values in the posterior more unlikely for smaller alphas?)
tic()
post_output <- alpha_inla(lis_Phi_, form, exposure, alphas, hist_plot=TRUE)
toc()

# hmmm the loglik looks funny... whyyyyyy
plot(alphas, unlist(post_output$marloglik), type="l")

sapply(post_output$summary.fitted.value, function(x) plot_pred(x[,1], NULL, PBC, compare=TRUE, count))

        # ====================================================================== #
        #                   Evaluate the result with some metrics                #
        # ====================================================================== #

# PLOT:: spatial prediction with the ridge glm selected lengthscale
# plot_pred(exp(stan_mod$linear.predictors))
plot_pred(colMeans(yrep_), round(best_alpha, 2), PBC)
bias(yrep_)
rmse(yrep_)
cp95(yrep_)

for(i in 1:length(alphas)){
        plot(plot_pred(post_output$summary.fitted.value[[i]][,1], alphas[i], PBC
                       # , compare=FALSE
                       ))
}
