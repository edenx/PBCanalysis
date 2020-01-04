library(sf)
library(raster)
library(tidyverse)
library(tictoc)
library(ggplot2)
library(dplyr)
library(reshape2)

library(randtoolbox)
library(glmnet)
library(rstan)
library(rstanarm)
library(INLA)

# set working directory to the main repo
setwd("~/Documents/PBCanalysis")
source("ArealKernel/Functions/WrapperFunc.R")
source("ArealKernel/Functions/ModelsFunc.R")

# load data
load("ArealKernel/SpatialData/PBCshp.Rdata")
load("ArealKernel/SpatialData/pop_nut.Rdata")

# ------------------------------------ Kernel Preprocessing -------------------------------------

sf_dat <- st_as_sf(PBCshp) 
# PBC <- sf_dat
df_dat <- sf_dat %>% st_set_geometry(NULL)
count <- df_dat$X
regcovar <- select(df_dat, Income, Crime, Environment, Employment, 
                   Barriers, propmale, Education)

# precompute the RFF (Phi matrix; basis) for each lengthscale
alphas <- seq(0.2, 1, length.out = 30)

tic()
all.lis_Phi <- precomp_ker(PBCshp, pop_nut, alphas=alphas, RFF=TRUE, m=100,
                       if_response=TRUE, response=count,
                        if_regcovar=TRUE, regcovar=regcovar)
toc()

lis_Phi.null <- all.lis_Phi[[1]]
lis_Phi.resp <- all.lis_Phi[[2]]
lis_Pho.covar <- all.lis_Phi[[3]]

# Full covariance matrix without approximation
# it takes a long time to compute the full kernel, but only needs to compute once
# if the CBCV directly set the left out region's entries as 0
tic()
all.lis_Ker <- precomp_ker(PBCshp, pop_nut, alphas=alphas, 
                           if_response=TRUE, response=count,
                           if_regcovar=TRUE, regcovar=regcovar)
toc()

# --------------------------------------- GLM with RFF -----------------------------------------
# Stan Fit with no covariates ------------------------------------------------------------------

options(mc.cores = parallel::detectCores())
SEED <- 727282

glm.stan_lis <- list()
glm.stan_pred_lis <- list()

tic.clearlog()
tic("model fitting")
for(i in 1:length(lis_Phi.null)){
        dat <- as.data.frame(lis_Phi.null[[i]]) %>%
                mutate(count=count
                       , pop=scale(df_dat$pop, center=FALSE)
                )
        
        tic(paste0("lengthscale=", alphas[i]))
        glm.stan_lis[[i]] <- stan_glm(count ~ . -pop,
                                  offset=pop,
                                  data=dat, family=poisson,
                                  # prior=normal(0, sqrt(2)), prior_intercept=normal(0,5),
                                  control = list(max_treedepth = 20),
                                  # chains=5, 
                                  # thin=5,
                                  seed=SEED, verbose=TRUE)
        toc(log=TRUE)
        
        glm.stan_pred_lis[[i]] <- posterior_predict(glm.stan_lis[[i]])
        
}
toc()

# time taken for model fitting
lis_elap <- tic.log(format=FALSE)
elap <- sapply(lis_elap, function(x) x$toc - x$tic)
plot(alphas, elap, type="l", xlab="Lengthscale", ylab="time elapsed/sec", 
     main="time taken to sample")

# inspect the model fitting for one lengthscale
shinystan::launch_shinystan(glm.stan_lis[[8]])
pred_glm_stan <- glm.stan_pred_lis[[8]]
plot_pred(pred_glm_stan$mean, ls=alphas[[8]], sf_dat, compare=TRUE, count=count)

# save the result to local file
saveRDS(glm.stan_lis, "~/Documents/StanResult/stan_mod_rff.RData") #0.09285714; 0.11357143
saveRDS(glm.stan_pred_lis, "~/Documents/StanResult/stan_pred_rff.RData")

# INLA Fit with no covariates ------------------------------------------------------------------

lmod <- lm(count ~ . , lis_Phi.null[[1]])
form <- formula(lmod)
# flat exposure for now
exposure <- rep_along(df_dat$pop, 1)

tic()
glm.inla_lis <- alpha_inla(lis_Phi.null, form, exposure, alphas, hist_plot=FALSE)
toc()

# hmmm the loglik looks funny... whyyyyyy
plot(alphas, unlist(glm.inla_lis$marloglik), type="l")
# sapply(glm.inla_lis$summary.fitted.value, function(x) plot_pred(x[,1], NULL, sf_dat, compare=TRUE, count))

# inspect the model fitting for one lengthscale
pred_glm_inla <- glm.inla_lis[[8]]
plot_pred(pred_glm_inla$summary.fitted.value[,1], alphas[8], sf_dat, compare=TRUE, count)

# ------------------------------------- MVN with Full Kernel: INLA ---------------------------------
# create aggregated Kernel for a specific lengthscale ls (full)
ls <- alphas[8]
regker <- create_ker(ls, PBCshp, pop_den, plot=TRUE)
inv_regker <- Matrix::solve(regker)

df_dat.bym2 <- df_dat %>% mutate(ID=1:nrow(.))
# flat exposure for now
exposure <- rep_along(df_dat$pop, 1)

form.covar <- as.formula(paste("X", paste(c("Income", "Crime", "Environment", 
                                          "Employment", "Barriers", "propmale", 
                                          "Education"), collapse=" + "),
                             sep=" ~ "))
form.null <- X ~ 1

# INLA fit with one lengthscale ----------------------------------------------------------------
tic()
bym3.inlafit <- inla(update(form.covar, . ~. +
                               f(ID, model="bym2", graph=inv_regker, 
                                 scale.model=TRUE, constr=TRUE)), 
                data=PBC_df.bym, 
                family="poisson", E=exposure, 
                control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE),
                control.predictor=list(compute = TRUE)
)
toc()
plot_pred(bym3.inlafit$summary.fitted.values[,1], ls, PBC, compare=TRUE, count)


#-------------------------------------- MVN with RFF: Stan -------------------------------------

# Stan fit with one lengthscale -----------------------------------------------------------------
options(mc.cores = parallel::detectCores())

N <- as.integer(nrow(PBC))
m <- 100
l <- as.integer(2*m)
count <- as.integer(count)
mat.desgn <- lis_Phi.null[[8]]
# flat exposure for now
exposure <- rep_along(df_dat$pop, 1)

tic()
bym3.stan <- stan_model("ArealKernel/SpatialRE/bym3.stan")
bym3.stanfit <- sampling(bym3.stan, data=list(N=N,y=count,E=exposure,l=l,L=mat.desgn), 
                     control = list(adapt_delta = 0.97), 
                     chains=3, warmup=6000, iter=8000, save_warmup=FALSE);
toc()
# 1651.485 sec elapsed

# inspect the model fitting
shinystan::launch_shinystan(bym3.stanfit)
print(bym3.stanfit, pars=c("lp__", "beta0", "rho", "sigma", "log_precision", 
                       "logit_rho", "mu[5]", "phi[10]", "theta[5]"), 
      probs=c(0.025, 0.5, 0.975))
pred_bym3_stan <- stan_interval(bym3.stanfit, "mu")
plot_pred(pred_bym3_stan$mean, alphas[8], sf_dat, compare=TRUE, count=count)

# save the fit
saveRDS(bym3.stanfit, file="~/Documents/StanResult/bym3_stanfit.Rdata")

# 
# # ----------------------------------------------------------------------------------------------
# #                                       BYM2 Model
# # ----------------------------------------------------------------------------------------------
# 
# dat.adj <- poly2nb(sf_dat)
# W.dat <- nb2mat(dat.adj, style="B")
# # W.dat.rs <- nb2mat(dat.adj, style="W")
# 
# df_dat.bym2 <- df_dat %>% mutate(ID=1:nrow(.))
# # flat exposure for now
# exposure <- rep_along(df_dat$pop, 1)
# 
# form.covar <- as.formula(paste("X", paste(c("Income", "Crime", "Environment", 
#                                             "Employment", "Barriers", "propmale", 
#                                             "Education"), collapse=" + "),
#                                sep=" ~ "))
# form.null <- X ~ 1
# 
# # INLA -----------------------------------------------------------------------------------------
# tic()
# # ICAR component as spatial effect
# bym2.inlafit <- inla(update(form.covar, . ~. +
#                                f(ID, model="bym2", graph=W.dat, 
#                                  scale.model=TRUE, constr=TRUE)), 
#                 data=df_dat.bym2, 
#                 family="poisson", E=exposure, 
#                 control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE),
#                 control.predictor=list(compute = TRUE)
# )
# toc()
# plot_pred(bym2.inlafit$summary.fitted.values[,1], "NULL", sf_dat, compare=TRUE, count)
# 
# # Stan -----------------------------------------------------------------------------------------
