library(SDALGCP)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(reshape2)

# Bias
bias <- function(pred) sapply(pred, function(x) mean(colMeans(x-count)))

# RMSE
rmse <- function(pred) sapply(pred, function(x) sqrt(mean(colMeans((x-count)^2))))

# Credible Interval 95%
cp95 <- function(pred){
        ci95 <- lapply(pred, function(x) apply(x, 2, quantile, c(0.975, 0.025)))
        sapply(ci95, function(x) mean(x[1,]>=count & x[2, ]<=count))
}

# prediction: spatial plot
plot_pred <- function(pred, compare=TRUE, count="X"){
        
        if(compare){
                sf_pred <- PBC[, c("geometry", paste0(count))] %>% 
                        mutate(pred=pred, normpop=lis_wcentr) %>% 
                        rename(count=paste0(count))
                
                pop_ <- ggplot() + 
                        geom_sf(data=sf_pred["count"], aes(fill=count), lwd=0.05) + 
                        scale_fill_distiller(palette="BuPu", direction=1, limits=c(0,8)) +
                        theme_linedraw() +
                        labs(x="Northing", y="Easting")
                
                pred_ <- ggplot() + 
                        geom_sf(data=sf_pred["pred"], aes(fill=pred), lwd=0.05) + 
                        scale_fill_distiller(palette="BuPu", direction=1, limits=c(0,8)) +
                        theme_linedraw() +
                        labs(x=paste0("Northing"), y="Easting")
                
                gridExtra::grid.arrange(pop_, pred_, nrow=2)
        }else{
                sf_pred <- PBC[, c("geometry")] %>% 
                        mutate(pred=pred, normpop=lis_wcentr)
                
                pred_ <- ggplot() + 
                        geom_sf(data=sf_pred["rsurf"], aes(fill=pred), lwd=0.05) + 
                        scale_fill_distiller(palette="BuPu", direction=1, limits=c(0,8)) +
                        theme_linedraw() +
                        labs(x=paste0("Northing"), y="Easting")
                
                print(rsurf_)
                
        }
}

# For incidence prediction
stan_mod <- readRDS("stan_mod_rff.RData") #0.09285714; 0.11357143
stan_mod_2 <- readRDS("stan_mod_wo.RData")
stan_mod_3 <- readRDS("stan_mod_w.RData")

yrep <- readRDS("stan_pred_rff.RData")
yrep_2 <- readRDS("stan_pred_wo.RData")
yrep_3 <- readRDS("stan_pred_w.RData")

PBC <- st_as_sf(PBCshp)
count <- PBC$X

## risk surface
rsurf <- lapply(stan_mod, function(x) exp(x$linear.predictors))
sapply(rsurf, plot_rsurf, compare=FALSE)

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

plot_pred <- rbind(bias_pred, rmse_pred, cp95_pred)
levels(plot_pred$model) <- c("Spatial risk with offset", 
                             "Covariates with/o offset",
                             "Covariates with offset")

ggplot(plot_pred, aes(alphas, value)) + 
        geom_line(aes(linetype=model)) +
        facet_grid(metric ~ ., scales="free") +
        theme(legend.position="top")




