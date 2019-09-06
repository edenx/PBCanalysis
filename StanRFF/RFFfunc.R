# for every intersected cells with polygon, normalise the population density per polygon
check.0 <- function(x) {
        w_val <- x[,"value"]*x[,"weight"]
        if(any(x[,"value"] != 0)){
                cbind(x, W=x[,"value"]/sum(w_val))}
        else{cbind(x, W=x[,"value"])}
}

rff.region <- function(wcentr, Omega, m, alpha){
        centr <- as.matrix(wcentr[,c("x", "y")])
        # population density empirical distribution
        w <- as.matrix(wcentr[,c("W")])
        h <- as.matrix(wcentr[,c("weight")])
        
        # Projection - combine data with sample frequencies
        proj <- centr %*% t(Omega) 
        
        # Fourier feature for a given area
        # print(w*h)
        phi <- sqrt(1/m) * colSums(t(w*h) %*% cbind(cos(proj/alpha), sin(proj/alpha)))
        
        return(phi)
}

# using the inversion -- quasi monte carlo integration;
# Matern 5-2 kernel -- student 5-2 sdf
sim_rff <- function(lis_region, m=100, nu=5, alpha, plot=FALSE){
        Omega <- qt(halton(m,2), nu)
        Phi <- t(sapply(lis_region, rff.region, Omega, m, alpha))
        
        Kernel <- Phi %*% t(Phi)
        if(plot){
                print(lattice::levelplot(Kernel))
        }
        
        return(Phi)
}

# find the best lengthscale with cv over ridge glm
alpha_cv <- function(lis_rff, offset, tune_param=1, hist_plot=TRUE){
        
        lis_cvfit <- c()
        lis_f <- list()
        lis_fit <- list()
        
        for(i in 1:length(lis_rff)){
                Phi <- lis_rff[[i]]
                print(i)
                
                # fit a regularised glm with log link
                lambdas <- 10^seq(5, -5, length.out=100)
                cv_fit <- cv.glmnet(Phi, count, family="poisson", 
                                    offset=log_pop,
                                    alpha=1, lambda=lambdas)
                
                opt_lambda <- cv_fit$lambda.min
                lis_cvfit <- c(lis_cvfit, cv_fit$cvm[which(cv_fit$lambda == opt_lambda)])
                
                cat("The optimal lambda is ", opt_lambda)
                cat("\n")
                
                fit <- glmnet(Phi, count, family="poisson", 
                              offset=log_pop,
                              alpha = 0.6, lambda=opt_lambda)
                lis_fit[[i]] <- fit
                
                lis_f[[i]] <- predict(fit, s=opt_lambda, newx=Phi, 
                                      newoffset=log_pop,
                                      type="response"
                )
                
                if(hist_plot){
                        hist(lis_f[[i]])
                }
                
        }
        
        plot(alphas, lis_cvfit, type="l", xlab="lengthscale", 
             main="Cross Validation Error over Lengthscales")
        
        min_index <- which.min(lis_cvfit)
        # the best lengthscale with ridge glm
        best_alpha <- alphas[min_index]
        best_pred <- lis_f[[min_index]]
        
        out <- list(min_index=min_index, best_ls=best_alpha, 
                    best_pred=best_pred, lis_cvfit=lis_cvfit)
        return(out)
}

# Bias
bias <- function(pred) mean(colMeans(pred-count))

# RMSE
rmse <- function(pred) sqrt(mean(colMeans((pred-count)^2)))

# Credible Interval 95%
cp95 <- function(pred){
        ci95 <- apply(pred, 2, quantile, c(0.975, 0.025))
        mean(ci95[1,]>=count & ci95[2, ]<=count)
}

# prediction: spatial plot
plot_pred <- function(pred, ls, compare=TRUE, count="X"){
        
        if(compare){
                sf_pred <- PBC[, c("geometry", paste0(count))] %>% 
                        mutate(pred=pred, normpop=lis_wcentr) %>% 
                        rename(count=paste0(count))
                
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
                sf_pred <- PBC[, c("geometry")] %>% 
                        mutate(pred=pred, normpop=lis_wcentr)
                
                pred_ <- ggplot() + 
                        geom_sf(data=sf_pred["pred"], aes(fill=pred), lwd=0.05) + 
                        scale_fill_distiller(palette="BuPu", direction=1) +
                        theme_linedraw() +
                        labs(x=paste0("Northing/ lengthscale=", ls), y="Easting")
                print(pred_)
        }
}