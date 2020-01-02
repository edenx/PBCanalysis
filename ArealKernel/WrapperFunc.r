source("RegKernelFunc.r")

## creat regional aggregated kernel --------------------------------------------------------------

# alpha, nu: real, hyperparameter of kernel (Matern)
# dat: sp or sf object that contains spatial polygons dataframe
# pop.rast: raster object that overlaps with dat
# RFF: if Random Fourier features approximation should be used
# m: if RFF is true, the number of random features used
# plot: if the heatmap of produced kernel should be plotted 

# NOTE: 
# 1. RFF approximation would result in non-psd kernel
# 2. The RFF option returns Cholesky decomposition of the cov matrix if plot=FALSE
#       otherwise the full covariance is returned
# 3. Exposure calculation is not implemented; is the comment valid?
create_ker <- function(alpha, dat, pop.rast, nu=5, RFF=FALSE, m=NULL, plot=FALSE){
        pop_den_ <- raster::intersect(pop.rast, dat) %>% replace_na(0)
        sf_dat <- st_as_sf(dat)

        # value: population density per cell
        # weight: the proportion of cell intersected with the polygon
        # W: weighted population density per grid intersection with polygon
        sf_poly <- raster::extract(pop_den_, sf_dat, cellnumbers=TRUE, small=TRUE,
                                    weights=TRUE, normalizeWeights=FALSE)
        lis_centr <- sapply(sf_poly, function(x) cbind(x, scale(coordinates(pop_den_))[x[,1],]))
        lis_wcentr <- sapply(lis_centr, check.0)
        
        # # find the exposure risk
        # pop_den_reg <- sapply(lis_centr, function(x) sum(x[,"value"] * x[,"weight"]))
        # exposure <- df_dat$pop/pop_den_reg
        
        # find the aggregated kernel
        if(RFF){
                regker <- sim.rff(lis_wcentr, alpha=alpha, m=m, nu=nu, plot=plot)
        }else{
                regker <- sim.regker(lis_wcentr, alpha=alpha, nu=nu, plot=plot)
        }

        return(regker)
}

## precompute the cov matrix for all lengthscales -----------------------------------------------

# alphas: real vector, lengthscale of the kernel (Matern)
# response: real vector; response value of spatial model
# regcov: dataframe; spatial covariates same level as in dat

# NOTE: "plot" option is supressed

precomp_ker <- function(dat, pop.rast, nu=5, alphas=seq(0.2, 1, length.out = 30),
                        RFF=FALSE, m=NULL, plot=FALSE,
                        if_response=FALSE, response=NULL,
                        if_regcovar=FALSE, regcovar=NULL
                        ){
        
        lis_Phi.null <- sapply(alphas, create_ker, dat, pop.rast, nu=5,
               RFF=RFF, m=m, plot=FALSE)
        lis_Phi <- list()
        
        if(if_regcovar && if_response){
                for(i in alphas){
                        lis_Phi[[i]] <- lis_Phi.null[[i]] %>% as.data.frame() %>%
                        dplyr::mutate(resp=response) %>%
                        cbind(regcovar)
                }
        }else if(if_response){
                for(i in alphas){
                        lis_Phi[[i]] <- lis_Phi.null[[i]] %>% as.data.frame() %>%
                        dplyr::mutate(resp=response)
                }
        }else if(if_regcovar){
                for(i in alphas){
                        lis_Phi_[[i]] <- lis_Phi.null[[i]] %>% 
                        cbind(regcovar)
                }
        }
        
        if(if_response || if_regcovar){
                return(lis_Phi)
        }else
                return(lis_Phi.null)
}

#-------------------------------- Prediction: spatial plot ----------------------------------------
plot_pred <- function(pred, ls, sf_geom, compare=FALSE, count=NULL){
        
        if(compare){
                sf_pred <- sf_geom[, "geometry"] %>% 
                        mutate(pred=pred, count=count)
                
                pop_ <- ggplot() + 
                        geom_sf(data=sf_pred["count"], aes(fill=count), lwd=0.05) + 
                        scale_fill_distiller(palette="BuPu", direction=1) +
                        theme_linedraw() +
                        labs(x="Northing", y="Easting")
                
                pred_ <- ggplot() + 
                        geom_sf(data=sf_pred["pred"], aes(fill=pred), lwd=0.05) + 
                        scale_fill_distiller(palette="BuPu", direction=1) +
                        theme_linedraw() +
                        labs(x=paste0("Northing/ lengthscale=", ls), y="Easting")
                
                gridExtra::grid.arrange(pop_, pred_, nrow=2)
        }else{
                sf_pred <- sf_geom[, c("geometry")] %>% 
                        mutate(pred=pred)
                
                pred_ <- ggplot() + 
                        geom_sf(data=sf_pred["pred"], aes(fill=pred), lwd=0.05) + 
                        scale_fill_distiller(palette="BuPu", direction=1) +
                        theme_linedraw() +
                        labs(x=paste0("Northing/ lengthscale=", ls), y="Easting")
                print(pred_)
        }
}
