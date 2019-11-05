# for every intersected cells with polygon, normalise the population density per polygon
check.0 <- function(x) {
        w_val <- x[,"value"]*x[,"weight"]
        if(any(x[,"value"] != 0)){
                cbind(x, W=x[,"value"]/sum(w_val))}
        else{cbind(x, W=x[,"value"])}
}

# distance matrix between cells within two regions A, B
vectorized_pdist <- function(A, B){
        an <- apply(A, 1, function(rvec) crossprod(rvec,rvec))
        bn <- apply(B, 1, function(rvec) crossprod(rvec,rvec))
        
        m <- nrow(A)
        n <- nrow(B)
        
        tmp <- matrix(rep(an, n), nrow=m) 
        tmp <- tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
        sqrt(tmp - 2*tcrossprod(A,B))
}

# alpha lengthscale (maybe a bit confusing, change afterwards)
# nu any in c(5, 3, 1, "inf")
Matern.ker <- function(coord1, coord2, alpha, nu){
        # dim: nrow(coord1) x nrow(coord2)
        d <- vectorized_pdist(coord1, coord2)
        
        if(!nu %in% c(1, 3, 5, "inf")){
                errorCondition("choose appropriate nu in k52, k32, k12, inf")
        }else if(nu==5){
                (1 + sqrt(5)*d/alpha + 5*d^2/(3*alpha^2)) * exp(-sqrt(5)*d/alpha)
        }else if(nu==3){
                (1 + sqrt(3)*u/alpha) * exp(-sqrt(3)*d/alpha)
        }else if(nu==1){
                exp(-d/alpha)
        }else{
                exp(-d^2/alpha^2)
        }
}

Region.ker <- function(wcentr1, wcentr2, alpha, nu=5){
        centr1 <- as.matrix(wcentr1[,c("x", "y")])
        centr2 <- as.matrix(wcentr2[,c("x", "y")])
        
        # population density empirical distribution
        weight1 <-  as.matrix(wcentr1[,c("W")]) * as.matrix(wcentr1[,c("weight")])
        weight2 <-  as.matrix(wcentr2[,c("W")]) * as.matrix(wcentr2[,c("weight")])
        
        base.ker <- Matern.ker(centr1, centr2, alpha, nu)
        
        # weighted average over region with normalised population density
        t(weight1) %*% base.ker %*% weight2
}


sim.regker <- function(lis_region, alpha, nu=5, plot=FALSE){
        # parallel the outer for loop
        require(parallel)
        no_cores <- detectCores() - 1
        cl <- makeCluster(no_cores, type="FORK")
        
        Kernel <- parSapply(cl, lis_region, 
               function(region) parSapply(lis_region, Region.ker, region, alpha, nu))
        stopCluster(cl)
        
        if(plot){
                print(lattice::levelplot(Kernel))
        }
        
        return(Kernel)
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
        
        return(as.matrix(Phi))
}

# find the best lengthscale with cv over ridge glm
alpha_ridge <- function(lis_rff, E, tune_param=1, hist_plot=TRUE){
        
        lis_cvfit <- c()
        lis_f <- list()
        lis_fit <- list()

        for(i in 1:length(lis_rff)){
                Phi <- lis_rff[[i]]
                print(i)
                
                # fit a regularised glm with log link
                lambdas <- 10^seq(5, -5, length.out=100)
                cv_fit <- cv.glmnet(Phi, count, family="poisson",
                                    offset=E,
                                    alpha=tune_param, lambda=lambdas)

                opt_lambda <- cv_fit$lambda.min
                lis_cvfit <- c(lis_cvfit, cv_fit$cvm[which(cv_fit$lambda == opt_lambda)])

                cat("The optimal lambda is ", opt_lambda)
                cat("\n")
                
                fit <- glmnet(Phi, count, family="poisson", 
                              offset=E,
                              alpha=tune_param, lambda=opt_lambda)
                
                lis_fit[[i]] <- fit
                
                lis_f[[i]] <- predict(fit, s=opt_lambda, newx=Phi, 
                                      newoffset=E,
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

# find the posterior with INLA
# lmod <- lm(count ~ . + I(pop)-pop, dat_)
# form<- formula(lmod)

alpha_inla <- function(lis_rff, form, E, alphas, hist_plot=TRUE){
        lis_f <- list()
        lis_fe <- list()
        lis_fv <- list()
        lis_mlik <- list()

        for(i in 1:length(lis_rff)){
                Phi <- lis_rff[[i]]
                
                # prior=normal(0, sqrt(2)), prior_intercept=normal(0,5)
                fit <- inla(form, family="poisson", data=Phi, 
                            control.fixed(prec.intercept=1/sqrt(2), prec=10),
                            control.predictor = list(compute=TRUE), E=E)
                lis_f[[i]] <- fit
                lis_fe[[i]] <- fit$summary.fixed
                lis_fv[[i]] <- fit$summary.fitted.values
                lis_mlik[[i]] <- fit$mlik[1]
                
                if(hist_plot){
                        # mean of predicted count
                        col1 <- adjustcolor(2, alpha.f = 0.3)
                        col2 <- adjustcolor(4, alpha.f = 0.3)
                        hist(lis_fv[[i]][,1], col=col1, 
                             main=paste0("Data vs. Predicted Mean Count with ls=", alphas[i]))
                        hist(Phi$count, col=col2, add=TRUE)
                        legend("topright", legend=c("Predicted Mean", "Count"),
                               fill=c(col1, col2), cex=0.8)
                }
        }
        out <- list(fit=lis_f, summary.fixed.effect=lis_fe, 
                    summary.fitted.value=lis_fv, marloglik=lis_mlik)
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

# prediction: histogram
plot_hist <- function(pred, alpha){
        col1 <- adjustcolor(2, alpha.f = 0.3)
        col2 <- adjustcolor(4, alpha.f = 0.3)
        hist(pred, col=col1, 
             main=paste0("Data vs. Predicted Mean Count with ls=", alpha))
        hist(PBC_df$X, col=col2, add=TRUE)
        legend("topright", legend=c("Predicted Mean", "Count"),
               fill=c(col1, col2), cex=0.8)
}

# prediction: spatial plot
plot_pred <- function(pred, ls, compare=FALSE, count=NULL){
        
        if(compare){
                sf_pred <- PBC[, "geometry"] %>% 
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
                sf_pred <- PBC[, c("geometry")] %>% 
                        mutate(pred=pred)
                
                pred_ <- ggplot() + 
                        geom_sf(data=sf_pred["pred"], aes(fill=pred), lwd=0.05) + 
                        scale_fill_distiller(palette="BuPu", direction=1) +
                        theme_linedraw() +
                        labs(x=paste0("Northing/ lengthscale=", ls), y="Easting")
                print(pred_)
        }
}
