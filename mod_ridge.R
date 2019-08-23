#### Run 'RFFregion.R' first to get Fourier Features

library(tidyverse)
library(tictoc)
library(glmnet)

# find the best lengthscale
lis_cvfit <- c()
lis_f <- list()
lis_fit <- list()

for(i in 1:length(lis_Phi)){
        Phi_ <- lis_Phi_[[i]]
        print(i)
        
        # fit a regularised glm with log link
        lambdas <- 10^seq(5, -5, length.out=100)
        cv_fit <- cv.glmnet(Phi_, count, family="poisson", 
                            # offset=log_pop,
                            alpha=0.6, lambda=lambdas)
        
        opt_lambda <- cv_fit$lambda.min
        lis_cvfit <- c(lis_cvfit, cv_fit$cvm[which(cv_fit$lambda == opt_lambda)])
        
        cat("The optimal lambda is ", opt_lambda)
        cat("\n")
        
        fit <- glmnet(Phi_, count, family="poisson", 
                      # offset=log_pop,
                      alpha = 0.6, lambda=opt_lambda)
        lis_fit[[i]] <- fit
        
        lis_f[[i]] <- predict(fit, s=opt_lambda, newx=Phi_, 
                              # offset=log_pop,
                              type="response"
        )
        
        hist(lis_f[[i]])
}

plot(alphas, lis_cvfit, type="line")

min_error <- which.min(lis_cvfit)
best_alpha <- alphas[min_error]
best_pred <- lis_f[[min_error]]
hist(best_pred, main=paste0("Best Lengthscale is ", best_alpha))
