#### Run 'RFFregion.R' first to get Fourier Features

library(tidyverse)
library(tictoc)
library(rstanarm)
library(ggplot2)

options(mc.cores = parallel::detectCores())
SEED <- 727282

        # ====================================================================== #
        #               HMC: Sample from Posterior distribution                  #
        # ====================================================================== #

stan_mod_3 <- list()
yrep_3 <- list()

tic.clearlog()
tic("model fitting")
for(i in 1:length(lis_Phi_)){
        dat <- as.data.frame(lis_Phi_[[i]]) %>%
                mutate(count=count
                       , log_pop=scale(log(PBC$pop), center=FALSE)
                )
        
        print(paste0("alpha=", alphas[i]))
        # find the starting position for the HMC
        # start_mod <- glm(count ~ . - log_pop, offset=log_pop,
        #                  family="poisson", data=dat)
        # start_beta <-coef(start_mod)
        # start_sigma2 <- mean(start_mod$residuals^2)
        # plot(start_beta, type="l")
        
        tic(paste0("lengthscale=", alphas[i]))
        stan_mod_3[[i]] <- stan_glm(count ~ . -log_pop,
                                  offset=log_pop,
                                  data=dat, family=poisson,
                                  prior=normal(0, sqrt(2)), prior_intercept=normal(0,5),
                                  control = list(max_treedepth = 20),
                                  # chains=5, 
                                  # thin=5,
                                  seed=SEED, verbose=TRUE)
        toc(log=TRUE)
        
        yrep_3[[i]] <- posterior_predict(stan_mod_3[[i]])
        
}
toc()

lis_elap <- tic.log(format=FALSE)
elap <- sapply(lis_elap, function(x) x$toc - x$tic)
plot(alphas, elap, type="l", xlab="Lengthscale", ylab="time elapsed/sec", 
     main="time taken to sample")



# shinystan::launch_shinystan(stan_mod[[7]])

# use sample mean
for(i in 1:length(yrep_3)){
        pred <- yrep_3[[i]]
        hist(count, main="posterior mean vs. original count")
        hist(round(colMeans(pred)), add=TRUE, col="red", main=paste0("alpha=", alpha[i]))
        print(table(count))
        print(table(round(colMeans(pred))))
        
        
        sf_pred <- PBC[, c("geometry", "X")] %>% mutate(pred=colMeans(pred),
                                                        normpop=lis_wcentr) %>% rename(count=X)
        
        pop_ <- ggplot() + 
                geom_sf(data=sf_pred["count"], aes(fill=count), lwd=0.05) + 
                scale_fill_distiller(palette="BuPu", direction=1, limits=c(0,8)) +
                theme_linedraw() +
                labs(x="Northing", y="Easting")
        
        pop_pred <- ggplot() + 
                geom_sf(data=sf_pred["pred"], aes(fill=pred), lwd=0.05) + 
                scale_fill_distiller(palette="BuPu", direction=1, limits=c(0,8)) +
                theme_linedraw() +
                labs(x=paste0("Northing/lengthscale=", alphas[i]), y="Easting")
        
        
        gridExtra::grid.arrange(pop_, pop_pred, nrow=2)
}

        # ====================================================================== #
        #               find lengthscale: kfold cross validation                 #
        # ====================================================================== #
# tic("5fold cv")
# lis_cv <- list()
# for(i in 1:length(stan_mod)){
#         lis_cv[[i]] <- rstanarm::kfold(stan_mod[[i]], K=5)
# }
# 
# lis_elpd_kfold <- sapply(lis_cv, function(x) x$elpd_kfold)
# toc()

library(parallel)

# Calculate the number of cores
no_cores <- detectCores() - 2

# Initiate cluster
cl <- makeCluster(no_cores, type="FORK")

# 5-fold cross validation
tic("5fold cv")

lis_cv <- parLapply(cl, stan_mod_3, rstanarm::kfold, K=5)
stopCluster(cl)
lis_elpd_kfold <- sapply(lis_cv, function(x) x$elpd_kfold)

toc()




# the smallest positive elpd_kfold gives the smallest error
plot(alphas, abs(lis_elpd_kfold), type="l")

# preferred model
index <- which.max(lis_elpd_kfold)
best_alpha <- alphas[index] # 878.9474 #0.8478947
best_mod <- stan_mod[[index]]
best_pred <- as.data.frame(yrep[[index]])

hist(count, main="posterior mean vs. original count")
hist(round(colMeans(best_pred)), add=TRUE, col="red", main=paste0("alpha=", alphas[i]))
print(table(count))
print(table(round(colMeans(best_pred))))


sf_pred <- PBC[, c("geometry", "X")] %>% mutate(pred=colMeans(best_pred),
                                                normpop=lis_wcentr) %>% rename(count=X)

pop_ <- ggplot() + 
        geom_sf(data=sf_pred["count"], aes(fill=count), lwd=0.05) + 
        scale_fill_distiller(palette="BuPu", direction=1, limits=c(0,8)) +
        theme_linedraw() +
        labs(x="Northing", y="Easting")

pop_pred <- ggplot() + 
        geom_sf(data=sf_pred["pred"], aes(fill=pred), lwd=0.05) + 
        scale_fill_distiller(palette="BuPu", direction=1, limits=c(0,8)) +
        theme_linedraw() +
        labs(x="Northing", y="Easting")

gridExtra::grid.arrange(pop_, pop_pred, nrow=2)



# # Create mode function.
# getmode <- function(v) {
#         uniqv <- unique(v)
#         uniqv[which.max(tabulate(match(v, uniqv)))]
# }
# 
# # Calculate the mode using the user function.
# best_mode <- lapply(best_pred, getmode)
# hist(unlist(best_mode))
