# for the RFF implementation, 
#       we require weoghted sum for the features w.r.t. polygon area/population density weights
# =============================================================================================

#       disease per unit area for each polygon
weighted.sim_lgcp <- function(mat){
        return(round(sum(apply(mat, 1, prod))))
}

#       total disease count for each polygon
sum.sim_lgcp <- function(mat){
        return(round(sum(mat[,1])))
}

#       polygon and pixel level population density statistics
sptpolyOpt <- function(rast, sf_poly, plot=TRUE){
        
        if(class(rast) != "RasterLayer") stop("Only RasterLayer object is considered.")
        if(any(class(sf_poly) != c('sf', 'data.frame'))) stop("Only sf object is considered.")
        
        # calculate weighted polulation of each spatial polygon by area
        #       i.e. reweighted population density for each polygon instead of cells
        sf_poly_ <- raster::extract(rast, sf_poly, cellnumbers=TRUE, 
                                fun=sum, sp=TRUE)
        sf_poly_ <- st_as_sf(sf_poly_) %>% rename(ArSum=NZ) 
        
        # find the overlaps of polygons and points (sf objects)
        sf_pixel <- as(as(rast, "SpatialGridDataFrame"), "sf")
        ints_polypnt <- st_intersection(sf_pixel, sf_poly_)
        
        # get (log) normalised (per polygon) population density 
        norm_pop_den <- aggregate(ints_polypnt$NZ, by=list(ints_polypnt$LSOA04CD), FUN=mean) %>%
                mutate(logmeanNZ=log(x)) %>%
                rename(LSOA04CD=Group.1, meanNZ=x) %>%
                mutate(sumNZ=aggregate(ints_polypnt$NZ, by=list(ints_polypnt$LSOA04CD), FUN=sum)$x)
        
        ints_polypnt <- ints_polypnt %>% left_join(norm_pop_den, by="LSOA04CD")  %>%
                mutate(normNZ=.$NZ/.$sumNZ)
        
        
        # sf objects with Polygons
        poly_polypnt <- aggregate(ints_polypnt, sf_poly_, unique)
        
        
        if(plot){
                # some plots
                area_ <- ggplot() + 
                        geom_sf(data=poly_polypnt[,"ArwNZ"], aes(fill=ArwNZ), lwd=0.05) + 
                        scale_fill_distiller(palette="BuPu", direction=1) +
                        theme_linedraw() +
                        labs(x="Northing", y="Easting", title="Cell proportion weighted aggregated pop_den")
                
                mean_ <- ggplot() + 
                        geom_sf(data=poly_polypnt[,"meanNZ"], aes(fill=meanNZ), lwd=0.05) + 
                        scale_fill_distiller(palette="BuPu", direction=1) +
                        theme_linedraw() +
                        labs(x="Northing", y="Easting", title="Poly mean pop_den")
                
                plot_norm <- as_tibble(st_coordinates(ints_polypnt)) %>% mutate(den=ints_polypnt$normNZ) 
                norm_ <- ggplot(plot_norm) + geom_raster(aes(x=X, y=Y, fill=den), interpolate=TRUE) + 
                        scale_fill_distiller(palette="BuPu", direction=1) +
                        geom_sf(data=poly_polypnt["geometry"], fill=NA, lwd=0.1) +
                        theme_linedraw() +
                        labs(x="Northing", y="Easting", fill="popwNZ", 
                             title="Total population density weighted pop_den")
                
                plot_rast <- as_tibble(st_coordinates(ints_polypnt)) %>% mutate(den=ints_polypnt$NZ) 
                rast_ <- ggplot(plot_rast) + geom_raster(aes(x=X, y=Y, fill=den), interpolate=TRUE) + 
                        scale_fill_distiller(palette="BuPu", direction=1) +
                        geom_sf(data=sf_poly["geometry"], fill=NA, lwd=0.1) +
                        theme_linedraw() +
                        labs(x="Northing", y="Easting", fill="NZ", title="pop_den")
                
                
                gridExtra::grid.arrange(area_, mean_, norm_, rast_, nrow=2)
        }
        
        out <- list(sf.pixel=ints_polypnt, sf.poly=poly_polypnt)
        
        return(out)
}

