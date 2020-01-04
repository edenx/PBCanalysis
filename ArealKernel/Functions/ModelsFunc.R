
#-------------------------------------- Ridge Regression -----------------------------------
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

#----------------------------------------- INLA ---------------------------------------------
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

