library(sf)
library(raster)
library(tidyverse)
library(tictoc)
library(ggplot2)
library(dplyr)
library(reshape2)

source("WrapperFunc.r")

# Data preprocessing ------------------------------------------------------------------------------

## Get population density data from FB HRSL
# geotiff_file <- "https://data.humdata.org/dataset/8c2c0b1f-66af-4a8e-b30e-59ad2249ee24/resource/d83a3bad-b72a-4e4e-9be9-93b4c654ac0f/download/population_mwi_2018-10-01.zip"
# download.file(geotiff_file, "pop_malawi_geotiff.zip")
# unzip("pop_malawi_geotiff.zip")

# # Aggregation
# tif_name <- 'population_mwi_2018-10-01.tif' 
# pop_malawi <- raster(tif_name)
# res(pop_malawi)
# # from 30m by 30m aggregate to 3km by 3km
# pop_malawi.aggre <- aggregate(pop_malawi, fact=100)
# save(pop_malawi.aggre, file="pop_malawi.aggre.Rdata")

# get population density of malawi
load("pop_malawi.aggre.Rdata")
res(pop_malawi.aggre)
plot(pop_malawi.aggre)

# get spatial data of malawi
df <- readRDS("~/Documents/sae/data/prev_malawi_2015.rds") %>% 
        st_as_sf() %>%
        mutate(y = est * n_obs, # y: total number of cases (ignores survey design)
               l_prev = qlogis(est), # l_prev: log prevalence estimate
               l_prev_se = se / (est * (1 - est))) # l_prev_se: calculated by delta method

# Other models see "/sae/03_sae.r"

# INLA: bym3 ----------------------------------------------------------------------------------
# create aggregated Kernel
ls <- 0.4
regker <- create_ker(ls, df, pop_malawi.aggre, plot=TRUE)
inv_regker <- Matrix::solve(regker)

# # precompute the RFF for each lengthscale
# # can specify regional covariates
# temp <- precomp_ker(df, pop_malawi.aggre)

inla_df <- list(y = round(df$y), m = df$n_obs, id=1:nrow(df))
inla_form <- y ~ 1 + f(id, model = "bym2", graph = inv_regker, scale.model = TRUE, constr=TRUE)

tic()
inla.fit_bym3 <- inla(inla_form,
                      family = "binomial",
                      control.family = list(control.link = list(model = "logit")),
                      data = inla_df, 
                      Ntrials = m,
                      control.predictor = list(compute = TRUE),
                      control.compute = list(dic = TRUE))
toc()
plot_pred(inla.fit_bym3$summary.fitted.values[,1], ls, df, compare=TRUE, df$est)
summary(inla.fit_bym3)



