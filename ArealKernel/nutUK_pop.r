library(sf)
library(raster)
library(tidyverse)
library(tictoc)
library(ggplot2)
library(dplyr)
library(reshape2)

library(randtoolbox)
library(SDALGCP)
library(glmnet)
library(rstanarm)
library(INLA)

setwd("~/Documents/PBCanalysis/ArealKernel")
source("WrapperFunc.r")

sf_dat <- st_as_sf(PBCshp)
df_dat <- sf_dat %>% st_set_geometry(NULL)
count <- df_dat$X
regcovar <- select(df_dat, Income, Crime, Environment, Employment, 
                   Barriers, propmale, Education)

# create aggregated Kernel
ls <- 0.4
regker <- create_ker(ls, PBCshp, pop_den, plot=TRUE)
inv_regker <- Matrix::solve(regker)

# precompute the RFF for each lengthscale
alphas <- seq(0.2, 1, length.out = 30)

# it takes a long time to compute the full kernel, but only needs to compute once
# if the CBCV directly set the left out region's entries as 0
tic()
lis_Phi <- precomp_ker(PBCshp, pop_den, alphas=alphas, if_response=TRUE, response=count,
                        if_regcovar=TRUE, regcovar=regcovar)
toc()



