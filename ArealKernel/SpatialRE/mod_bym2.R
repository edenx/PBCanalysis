library(rgdal)
library(sp)
library(spdep)
library(SDALGCP)
library(sf)

source("RFFfunc.R")

pop_den_ <- raster::intersect(pop_den, PBCshp) %>% replace_na(0)
PBC <- st_as_sf(PBCshp)

PBC_df.bym <- PBC %>% st_set_geometry(NULL) %>% 
        mutate(ID=1:nrow(.))

PBC.adj <- poly2nb(PBC)
W.PBC <- nb2mat(PBC.adj, style="B")
W.PBC.rs <- nb2mat(PBC.adj, style="W")

PBC.form <- as.formula(paste("X", paste(c("Income", "Crime", "Environment", 
                                          "Employment", "Barriers", "propmale", 
                                          "Education"), collapse=" + "),
                             sep=" ~ "))
PBC.form_ <- X ~ 1

tic()
PBC.bym <- inla(update(PBC.form, . ~. +
                                  f(ID, model="bym2", graph=W.PBC, 
                                    scale.model=TRUE, constr=TRUE)), 
                data=PBC_df.bym, 
                family="poisson", E=exposure, 
                control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE),
                control.predictor=list(compute = TRUE)
                   )
toc()

plot_pred(PBC.bym$summary.fitted.values[,1], "NULL", PBC, compare=TRUE, count)


# # Fitting instead with `gemeric3` r.e.
# Q  <- inla.scale.model(diag(rowSums(M)) - M,
#                        constr=list(A = matrix(1, 1, nrow(M)), e=0))
# tic()
# PBC.gen3 <- inla(update(PBC.form, . ~. +
#                                f(ID, model="generic3", 
#                                  Cmatrix=list(Q+diag(nrow(W.PBC)), Q),
#                                  hyper=list(theta11=list(initial=log(1/30000))))), 
#                 data=PBC_df.bym, 
#                 family="poisson", E=exposure, 
#                 control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE),
#                 control.predictor=list(compute = TRUE)
# )
# toc()
# summary(PBC.gen3)
# 
# plot_pred(PBC.gen3$summary.fitted.values[,1], "NULL", PBC, compare=TRUE, count)
# 
