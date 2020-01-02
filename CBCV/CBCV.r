# Inversely weighted by the number of neighbours (to avoid too many blocks to be selected)?
# Any justification for that?
library(spdep)
library(SDALGCP)
library(sf)
library(tidyverse)
library(dplyr)
library(reshape2)
library(ggplot2)
source("RFFfunc.R")

pop_den_ <- raster::intersect(pop_den, PBCshp) %>% replace_na(0)
PBC <- st_as_sf(PBCshp)
PBC_df <- PBC %>% st_set_geometry(NULL)
# ---------------------------------------------------------------------------------------- #
#                               Extending LOO to CBLOO       
# ---------------------------------------------------------------------------------------- #
# Extending LOO to CBLOO
# If randomly drawn from set of regions, set index=NULL by default
# If deterministic (i.e. loop over entire set of regions), index must be taken
CBLOO <- function(dat, rand=TRUE, index=NULL){
        mat_neighbour <- poly2nb(dat)
        if(rand){
                samp <- sample(1:nrow(dat),1)
        }else{
                samp <- index
        }
        samp_neighbour <- mat_neighbour[[samp]]
        
        # replace the response for selected regions by NA
        dat_new <- dat 
        dat_new[c(samp, samp_neighbour), "X"] <- rep(NA, length(samp_neighbour)+1)
        
        return(dat_new)
}
blckrm.1 <- CBLOO(PBC)
plot_pred(blckrm.1$X, NULL, PBC)

# fit the model in INLA (this gives immediate prediction for the NAs)
# CV error for prediction of the left out block

# ---------------------------------------------------------------------------------------- #
#                       Extending CBLOO by leaving out x% regions       
# ---------------------------------------------------------------------------------------- #

# ---------------------------- sample without replacement ---------------------------------
# nabla is the buffer for selecting roughly the integer number of regions given by the proportion
CBLOO_perc <- function(x, dat, nabla=2){
        mat_neighbour <- poly2nb(dat)
        dat_rmgeom <- dat %>% st_set_geometry(NULL) %>% mutate(index=1:n)
        
        samp_neighbour <- c()
        while(length(samp_neighbour)<round(nrow(dat)*x)-nabla){
                # sample one region and its neighbours
                samp <- sample(dat_rmgeom[,"index"],1)
                if(!samp %in% samp_neighbour){
                        samp_neighbour <- c(samp_neighbour, samp, mat_neighbour[[samp]])
                }
        }
        # update the data
        dat_new <- dat
        dat_new[samp_neighbour, "X"] <- rep(NA, length(samp_neighbour))
        
        # return(length(samp_neighbour))
        return(dat_new)
}
blckrm.05 <- CBLOO_perc(0.05, PBC)
plot_pred(blckrm.05$X, NULL, PBC)

# -------------------------------- sample with 1 draw ---------------------------------
# lambda is the number of samples to be drawn for estimating the number of neibours
# how do we make sure that the proportion is roughly what we want?
        # 
CBLOO_perc_1 <- function(x, dat, lambda=10){
        mat_neighbour <- poly2nb(dat)
        dat_rmgeom <- dat %>% st_set_geometry(NULL) %>% mutate(index=1:n)
        
        # find the number of neighbours for each region
        freq_neighbour <- sapply(mat_neighbour, length)
        # find the number of regions to directly sample from
        # take mean of random draw of lambda as the estimate of the number of neighbours
        samp_reg <- nrow(dat)*x/round(mean(sample(freq_neighbour,lambda))+1)
        
        # sample from the regions and take along with the neighbours
        samp <- sample(dat_rmgeom[,"index"],round(samp_reg))
        samp_neighbour <- c(samp, unlist(sapply(samp, function(x) mat_neighbour[[x]])))
        samp_neighbour <- unique(samp_neighbour)
        
        # update the data
        dat_new <- dat
        dat_new[samp_neighbour, "X"] <- rep(NA, length(samp_neighbour))
        
        # return(length(samp_neighbour))
        return(dat_new)
}

# visualise the new data
blckrm.05 <- CBLOO_perc_1(0.05, PBC, 10)
plot_pred(blckrm.05$X, NULL, PBC)

# # Can verify that the the mean of the number of left out reg is around x*n
# len_rm <- sapply(rep(0.05, 100), CBLOO_perc, PBC, 10)
# hist(len_rm)
# abline(v=mean(len_rm))


# ---------------------------------- Stratified scheme -----------------------------------
# Sample with 1 draw as above
CBLOO_perc_stra <- function(x, dat, lambda=10, stratified=FALSE, level=NULL){
        mat_neighbour <- poly2nb(dat)
        dat_rmgeom <- dat %>% st_set_geometry(NULL) %>% mutate(index=1:n)
        
        freq_neighbour <- sapply(mat_neighbour, length)
        samp_reg <- nrow(dat)*x/round(mean(sample(freq_neighbour,lambda))+1)
        
        if(stratified){
                # reconstruct the district name and the district code into two cols
                # "district name" "code"
                split_char <- strsplit(as.character(level), '[[:space:]]')
                lab_reg <- sapply(split_char, 
                                  function(x) c(paste0(x[-length(x)], collapse=""),x[length(x)]))
                lab_reg_df <- t(lab_reg) %>% as.data.frame(.) %>% 
                        rename("District"=V1, "Code"=V2) %>%
                        mutate("index"=1:nrow(.))
                
                # stratified sampling 
                # proportional to the number of regions within a district
                prob_grp <- lab_reg_df %>% group_by(District, add=TRUE) %>% 
                        summarise(freq=length(index))
                lab_reg_df_ <- merge(lab_reg_df, prob_grp, "District") %>% 
                        mutate(prob=freq/nrow(.))
                lab_reg_df_ <- lab_reg_df_[order(lab_reg_df_$index),]
                
                samp <- sample(lab_reg_df_$index, round(samp_reg), prob=lab_reg_df_$prob)
                
        }else{
                samp <- sample(dat_rmgeom$index, round(samp_reg))
        }
        samp_neighbour <- c(samp, unlist(sapply(samp, 
                                                function(x) mat_neighbour[[x]])))
        samp_neighbour <- unique(samp_neighbour)
        
        # update the data
        dat_new <- dat
        dat_new[samp_neighbour, "X"] <- rep(NA, length(samp_neighbour))
        
        return(dat_new)
}

# visualise the new data
blckrm.05 <- CBLOO_perc_stra(0.10, PBC, lambda=10, stratified=TRUE, level=PBC$LSOA04NM)
plot_pred(blckrm.05$X, NULL, PBC)
