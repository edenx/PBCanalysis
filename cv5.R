library(parallel)

# Calculate the number of cores
no_cores <- detectCores() - 2

# Initiate cluster
cl <- makeCluster(no_cores, type="FORK")

# 5-fold cross validation
tic("5fold cv")
lis_cv <- parLapply(cl, stan_mod_3, rstanarm::kfold, K=5)
stopCluster(cl)
lis_elpd_kfold <- sapply(lis_cv, function(x) x$elpd_kfold)
toc()


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
