library(SDALGCP)
library(sf)
library(tidyverse)
library(dplyr)
library(tictoc)
library(INLA)

source("RFFfunc.R")

pop_den_ <- raster::intersect(pop_den, PBCshp) %>% replace_na(0)
PBC <- st_as_sf(PBCshp)
PBC_df <- PBC %>% st_set_geometry(NULL)
exposure <- PBC_df$pop/pop_den_reg

# Can we incorporate bym2 framework into INLA for the derived spatial random effect f~MVN ???
# Can use Stan tho...
# worth a try...

## Naive implementation below: individual iid and spatial r.e., 
## with no common param controlling variance
## need to figure out the scaling and prior selection

inla_reg.meff <- function(Phi, mat_desn, form, Cmat, offset_){
        fit <- inla(
                update(form,. ~. + f(ID, model = "z", Z=Phi, Cmatrix=Cmat) 
                       + f(ID_c, model="iid")
                       ),
                data=mat_desn, family = "poisson", E=offset_,
                control.predictor = list(compute=TRUE))
        
        fv <- fit$summary.fitted.values
        mlik <- fit$mlik[1]
        
        out <- list(fv=fv, mlik=mlik, fit=fit)
        
        return(out)
}

# use mixed random effect for INLA
PBC_df.desn <- PBC_df %>% 
        mutate(ID=1:nrow(.), ID_c=1:nrow(.), count=.$X) %>%
        dplyr::select(ID, ID_c, count, Income, Crime, Environment, Employment, 
                      Barriers, propmale, Education)
# # offset
# pop <- scale(PBC_df$pop, center=FALSE)

# model formula
lmod <- lm(count ~ .-ID-ID_c, PBC_df.desn)
form <- formula(lmod)
Cmat <- diag(rep(1, ncol(lis_Phi.null[[1]])))
# try fitting with lengthscale=0.4
mat_desn <- PBC_df.desn
offset_ <- exposure
Phi <- lis_Phi.null[[which(abs(alphas-0.4)<0.01)]]
fit.4 <- inla(
        update(form,. ~. + f(ID, model = "z", Z=Phi, Cmatrix=Cmat) 
               + f(ID_c, model="iid")
        ),
        data=mat_desn, family = "poisson", E=offset_,
        control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE),
        control.predictor = list(compute=TRUE))
plot_pred(fit.4$summary.fitted.values[,1], 0.4, PBC, compare=TRUE, count)



# model fitting with multiple lengthscale
tic()
require(parallel)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, type="FORK")

lis_fit.inla <- parLapply(cl, lis_Phi.null, inla_reg.meff, PBC_df.desn, form, Cmat, exposure)

stopCluster(cl)
toc()

# # spatial plot with the prediction
# sapply(lis_fit.inla, function(x) plot_pred(x$fv[,1], NULL))

# find the marginal log lik of posterior
lis_mlik <- c()
for(i in 1:length(lis_fit.inla)){
        plot(lis_fit.inla[[i]]$fv[,1]-count)
        lis_mlik <- c(lis_mlik, lis_fit.inla[[i]]$mlik)
}
plot(alphas, lis_mlik, type="l")

# find the largest mlik
best_index <- which.max(lis_mlik)
best_alpha <- alphas[best_index]
# look at the spatial plot
plot_pred(lis_fit.inla[[best_index]]$fv[,1], best_alpha, compare=TRUE, count)


lapply(lis_fit.inla, function(x) table(round(x$fv[,1])))[[best_index]]
