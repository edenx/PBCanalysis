library(tidyverse)
library(tictoc)
library(rstanarm)
library(ggplot2)

setwd("~/Documents/PBCanalysis")
source("ArealKernel/Functions/WrapperFunc.R")

options(mc.cores = parallel::detectCores())
SEED <- 727282

stan_mod <- list()
yrep <- list()

tic.clearlog()
tic("model fitting")
for(i in 1:length(lis_Phi)){
        dat <- as.data.frame(lis_Phi[[i]]) %>%
                mutate(count=count
                       , log_pop=scale(log(PBC$pop), center=FALSE)
                )
        
        tic(paste0("lengthscale=", alphas[i]))
        stan_mod[[i]] <- stan_glm(count ~ . -log_pop,
                                  offset=log_pop,
                                  data=dat, family=poisson,
                                  prior=normal(0, sqrt(2)), prior_intercept=normal(0,5),
                                  control = list(max_treedepth = 20),
                                  # chains=5, 
                                  # thin=5,
                                  seed=SEED, verbose=TRUE)
        toc(log=TRUE)
        
        yrep[[i]] <- posterior_predict(stan_mod[[i]])
        
}
toc()

lis_elap <- tic.log(format=FALSE)
elap <- sapply(lis_elap, function(x) x$toc - x$tic)
plot(alphas, elap, type="l", xlab="Lengthscale", ylab="time elapsed/sec", 
     main="time taken to sample")

# save the result to local file
saveRDS(stan_mod, "~/Documents/StanResult/stan_mod_rff.RData") #0.09285714; 0.11357143
# saveRDS(stan_mod_2, "~/Documents/StanResult/stan_mod_wo.RData")
# saveRDS(stan_mod_3, "~/Documents/StanResult/stan_mod_w.RData")

saveRDS(yrep, "~/Documents/StanResult/stan_pred_rff.RData")
# saveRDS(yrep_2, "~/Documents/StanResult/stan_pred_wo.RData")
# saveRDS(yrep_3, "~/Documents/StanResult/stan_pred_w.RData")


# with the lengthscale selected by ridge glm
options(mc.cores = parallel::detectCores())
SEED <- 727282

dat <- as.data.frame(lis_Phi[[min_index]]) %>%
        mutate(count=count
               , log_pop=scale(log(PBC$pop), center=FALSE)
        )

tic(paste0("Model fitting with lengthscale=", alphas[min_index]))

stan_mod <- stan_glm(count ~ . -log_pop,
                          offset=log_pop,
                          data=dat, family=poisson,
                          prior=normal(0, sqrt(2)), prior_intercept=normal(0,5),
                          control = list(max_treedepth = 20),
                          # chains=5, 
                          # thin=5,
                          seed=SEED, verbose=TRUE)
toc()

yrep <- posterior_predict(stan_mod)


# when including the covariates
# with the lengthscale selected by ridge glm
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

