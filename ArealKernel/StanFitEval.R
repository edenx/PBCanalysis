library(ggplot2)
library(tidyverse)
library(dplyr)
library(reshape2)

setwd("~/Documents/PBCanalysis")
source("ArealKernel/Functions/WrapperFunc.R")
load("ArealKernel/SpatialData/PBCshp.Rdata")

# For incidence prediction
stan_mod <- readRDS("~/Documents/StanResult/stan_mod_rff.RData") #0.09285714; 0.11357143
stan_mod_2 <- readRDS("~/Documents/StanResult/stan_mod_wo.RData")
stan_mod_3 <- readRDS("~/Documents/StanResult/stan_mod_w.RData")

yrep <- readRDS("~/Documents/StanResult/stan_pred_rff.RData")
yrep_2 <- readRDS("~/Documents/StanResult/stan_pred_wo.RData")
yrep_3 <- readRDS("~/Documents/StanResult/stan_pred_w.RData")

PBC <- st_as_sf(PBCshp)
count <- PBC$X

## risk surface
rsurf <- lapply(stan_mod, function(x) exp(x$linear.predictors))
sapply(rsurf, plot_pred, compare=FALSE)

lis_pred <- list(yrep, yrep_2, yrep_3)

bias_pred <- as.data.frame(sapply(lis_pred, bias)) %>% 
        mutate(alphas=alphas) %>% 
        melt(id.vars="alphas", variable.name="Model") %>%
        mutate(metric=rep("Bias", nrow(.)))
        
rmse_pred <- as.data.frame(sapply(lis_pred, rmse)) %>%
        mutate(alphas=alphas) %>% 
        melt(id.vars="alphas", variable.name="Model") %>%
        mutate(metric=rep("RMSE", nrow(.)))

cp95_pred <- as.data.frame(sapply(lis_pred, cp95)) %>% 
        mutate(alphas=alphas) %>% 
        melt(id.vars="alphas", variable.name="Model") %>%
        mutate(metric=rep("CP95%", nrow(.)))

plot_pred_all <- rbind(bias_pred, rmse_pred, cp95_pred)
levels(plot_pred_all$model) <- c("Spatial risk with offset", 
                             "Covariates with/o offset",
                             "Covariates with offset")

ggplot(plot_pred_all, aes(alphas, value)) + 
        geom_line(aes(linetype=model)) +
        facet_grid(metric ~ ., scales="free") +
        theme(legend.position="top")


# PLOT:: spatial prediction with the ridge glm selected lengthscale
# plot_pred(exp(stan_mod$linear.predictors))
plot_pred(colMeans(yrep_), round(best_alpha, 2))
bias(yrep_)
rmse(yrep_)
cp95(yrep_)
