library(rgdal)
library(sp)
library(spdep)
library(SDALGCP)
library(sf)
library(INLA)

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
# ICAR component as spatial effect
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


tic()
# using the MVN with aggregated kernel as the spatial effect
PBC.ker <- inla(update(PBC.form, . ~. +
                                  f(ID, model="bym2", graph=inv_regker, 
                                    scale.model=TRUE, constr=TRUE)), 
                data=PBC_df.bym, 
                family="poisson", E=exposure, 
                control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE),
                control.predictor=list(compute = TRUE)
                   )
toc()
plot_pred(PBC.ker$summary.fitted.values[,1], "NULL", PBC, compare=TRUE, count)
