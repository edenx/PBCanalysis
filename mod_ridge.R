#### Run 'RFFregion.R' first to get Fourier Features
library(tidyverse)
library(tictoc)
library(glmnet)

source("RFFfunc.R")

log_pop <- log(PBC$pop)
# the best lengthscale with ridge glm
cv_output <- alpha_cv(lis_Phi_, tune_param=0.5, offset=log_pop)
min_index <- cv_output$min_index
best_alpha <- cv_output$best_ls
best_pred <- cv_output$best_pred

hist(best_pred, main=paste0("Best Lengthscale is ", best_alpha))
