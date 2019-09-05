## Not working!!!! WHY!!!
library(parallel)

# Calculate the number of cores
no_cores <- detectCores()

# Initiate cluster
# cl <- makeCluster(no_cores, type="FORK")

# lis_cv <- list()
# 5-fold cross validation
options(mc.cores=no_cores)
tic("5fold cv")
lis_cv <- lapply(stan_mod, rstanarm::kfold, K=5)

# for(i in 1:length(stan_mod_3)){
#         lis_cv[[i]] <- kfold(stan_mod_3[[i]], K=10, cores=8)
#         # lis_cv[[i]] <- rstanarm::loo(stan_mod_3[[i]], cores=no_cores)
#         paste0("current progress at", i)
# }

# stopCluster(cl)
lis_elpd_kfold <- sapply(lis_cv, function(x) x$elpd_kfold)
toc()

saveRDS(lis_cv, "stan_cv5.RData")

# the smallest positive elpd_kfold gives the smallest error
plot(alphas, abs(lis_elpd_kfold), type="l")

# preferred model
index <- which.max(lis_elpd_kfold)
best_alpha <- alphas[index] # 878.9474 #0.8478947
best_mod <- stan_mod[[index]]
best_pred <- as.data.frame(yrep[[index]])

hist(count, main="posterior mean vs. original count")
hist(round(colMeans(best_pred)), add=TRUE, col="red", main=paste0("alpha=", alphas[i]))
print(table(count))
print(table(round(colMeans(best_pred))))


sf_pred <- PBC[, c("geometry", "X")] %>% mutate(pred=colMeans(best_pred),
                                                normpop=lis_wcentr) %>% rename(count=X)

pop_ <- ggplot() + 
        geom_sf(data=sf_pred["count"], aes(fill=count), lwd=0.05) + 
        scale_fill_distiller(palette="BuPu", direction=1, limits=c(0,8)) +
        theme_linedraw() +
        labs(x="Northing", y="Easting")

pop_pred <- ggplot() + 
        geom_sf(data=sf_pred["pred"], aes(fill=pred), lwd=0.05) + 
        scale_fill_distiller(palette="BuPu", direction=1, limits=c(0,8)) +
        theme_linedraw() +
        labs(x="Northing", y="Easting")

gridExtra::grid.arrange(pop_, pop_pred, nrow=2)
