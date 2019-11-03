library(rgdal)
library(sp)
library(spdep)
library(SDALGCP)
library(sf)

source("RFFfunc.R")

pop_den_ <- raster::intersect(pop_den, PBCshp) %>% replace_na(0)
PBC <- st_as_sf(PBCshp)

PBC_df <- PBC %>% st_set_geometry(NULL) %>% 
        mutate(ID=1:nrow(.), pop=scale(.$pop, center=FALSE))

PBC.adj <- poly2nb(PBC)
W.PBC <- nb2mat(PBC.adj, style="B")
W.PBC.rs <- nb2mat(PBC.adj, style="W")

PBC.form <- as.formula(paste("X", paste(c("Income", "Crime", "Environment", 
                                          "Employment", "Barriers", "propmale", 
                                          "Education", "offset(pop)"), collapse=" + "),
                             sep=" ~ "))
PBC.form_ <- X ~ offset(pop)

PBC.bym <- inla(update(PBC.form, . ~. +
                                  f(ID, model = "bym2", graph = W.PBC)), 
                data=PBC_df, 
                family="poisson", 
                control.compute=list(dic = TRUE, waic = TRUE, cpo = TRUE),
                control.predictor=list(compute = TRUE)
                   )

plot_pred(PBC.bym$summary.fitted.values[,1], "NULL")
