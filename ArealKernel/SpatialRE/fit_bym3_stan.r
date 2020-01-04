library(rstan)
options(mc.cores = parallel::detectCores())

N <- as.integer(nrow(PBC))
m <- 100
l <- as.integer(2*m)
count <- as.integer(count)
mat.desgn <- lis_Phi.null[[8]]

tic()
bym3.stan <- stan_model("SpatialRE/bym3.stan")
bym3.fit <- sampling(bym3.stan, data=list(N=N,y=count,E=exposure,l=l,L=mat.desgn), 
                    control = list(adapt_delta = 0.97), 
                    chains=3, warmup=6000, iter=8000, save_warmup=FALSE);
toc()
# 1651.485 sec elapsed

# inspect the model fitting
shinystan::launch_shinystan(bym3.fit)
print(bym3.fit, pars=c("lp__", "beta0", "rho", "sigma", "log_precision", 
                                 "logit_rho", "mu[5]", "phi[10]", "theta[5]"), 
      probs=c(0.025, 0.5, 0.975))
pred_bym3 <- stan_interval(bym3.fit, "mu")
plot_pred(pred_bym3$mean, ls=NULL, PBC, compare=TRUE, count=count)

# save the fit
saveRDS(bym3.fit, file="~/Documents/StanResult/bym3_stanfit.Rdata")


