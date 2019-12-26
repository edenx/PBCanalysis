library(randtoolbox)
library(SDALGCP)
library(sf)
library(raster)
library(tidyverse)

source("RFFfunc.R")

pop_den_ <- raster::intersect(pop_den, PBCshp) %>% replace_na(0)
PBC <- st_as_sf(PBCshp)

sf_poly_ <- raster::extract(pop_den_, PBC, cellnumbers=TRUE, small=TRUE,
                            weights=TRUE, normalizeWeights=FALSE)
lis_centr <- sapply(sf_poly_, function(x) cbind(x, scale(coordinates(pop_den_))[x[,1],]))
lis_wcentr <- sapply(lis_centr, check.0)

# given lengthscale 0.6
test_phi <- sim_rff(lis_wcentr, alpha=0.6, plot=TRUE)
test_cov <- test_phi %*% t(test_phi)

mat_09 <- matrix(nrow=545, ncol=545)

# filter out for covariance greater than 0.8 
mat_09[which(test_cov > 0.8)] <- 1
mat_09[is.na(mat_09)] <- 0
lattice::levelplot(mat_09)

# diagonal approach: find the right diagonal sum (just as plotted with levelplot(.))
# should hit roughly zero values if encouters gaps between blocks
rot_mat <- t(apply(mat_09, 2, rev))

diag_mat <- split(rot_mat, col(rot_mat) - row(rot_mat))
temp_lis <- sapply(diag_mat, sum)

pnts_filtered <- rep_along(NA, temp_lis)
pnts_filtered[which(temp_lis==0)] <- temp_lis[which(temp_lis==0)]

# the trend of diagonal sum
plot(names(temp_lis), temp_lis, type="l")
points(names(temp_lis), pnts_filtered, col="red", pch=19, cex=0.2)

temp_index <- as.numeric(names(temp_lis[which(temp_lis==0)])) + 545
temp_round <- ceiling(1 + 0.5*(temp_index-1))
plot(diff(temp_round), type='l')   

sel_index <- temp_round[which(diff(temp_round) > 10)]
fin_index <- 545-c(1,temp_round[c(4,5,6,7,9,10,11,12,14,15)])
# this is the selected final indeces
fin_index <- sort(fin_index[-c(11,8)])
fin_index <- c(1, fin_index)

fin_index <- c(1,117,130,283,296,397,408,499,524,536,544)

mat_00 <- mat_09
for(i in fin_index){
        mat_00[0:i,i] <- 0.5
        mat_00[i,0:i] <- 0.5
}

png(paste0("CovBlock", ".png"), width=500)
lattice::levelplot(mat_00)
dev.off()

# # row approach
# lis_break <- c()
# mat_cur <- mat_09
# 
# find_block <- function(lis, prop){
#         if(mean(lis) < prop){
#                 return(i)
#         }
# }
# 
# mat_up <- mat_09
# mat_up[lower.tri(mat_up)] <- 0
# lis_index <- c()
# lis_mat <- list()
# n <- 1
# while(nrow(mat_up) > 2){
#         lis_prop <- apply(mat_up, 2, sum) / seq_along(mat_up[1,])
#         
#         dis_index <- which(lis_prop<0.5)[1]
#         mat_up <- mat_up[dis_index:nrow(mat_up), dis_index:nrow(mat_up)]
#         
#         lis_index <- c(lis_index, dis_index)
#         n <- n+1
#         lis_mat[[n]] <- mat_up
# }
# 
# length(lis_mat)
