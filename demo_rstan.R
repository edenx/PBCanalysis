library(rstanarm)
options(mc.cores = parallel::detectCores())
SEED <- 727282

###### try again with alpha=0.56
dat <- as.data.frame(lis_Phi_[[6]]) %>%
        mutate(count=count
               , log_pop=log(PBC$pop)
        )

# find the starting position for the HMC
start_mod <- glm(count ~ .-log_pop, 
                 offset=log_pop,
                 family="poisson", data=dat)
start_beta <- coef(start_mod)
start_sigma2 <- mean(start_mod$residuals^2)
# plot(start_beta, type="l")

tic("Stan")
stan_glm1 <- stan_glm(count ~ .-log_pop, 
                      offset=log_pop,
                      data=dat, family=poisson, 
                      prior=normal(start_beta[-1], 0.8), prior_intercept=normal(0,5),
                      control = list(max_treedepth = 20),
                      # chains=5, 
                      # thin=5,
                      seed=SEED, verbose=TRUE)
toc(log=TRUE)

yrep <- posterior_predict(stan_glm1)

library(shinystan)
shinystan::launch_shinystan(stan_glm1)

# predict

dat.yrep <- as.data.frame(yrep)
hist(count)
hist(round(colMeans(dat.yrep)), add=TRUE, col="red")
table(count)
table(round(colMeans(yrep)))

sf_pred <- PBC[, c("geometry", "X")] %>% mutate(pred=round(colMeans(yrep))) %>% rename(count=X)

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
