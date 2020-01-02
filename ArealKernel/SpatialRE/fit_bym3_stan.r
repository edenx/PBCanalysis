library(rstan)
options(mc.cores = parallel::detectCores())

N <- as.integer(nrow(PBC))
l <- as.integer(2*100)
count <- as.integer(count)
mat.desgn <- lis_Phi.null[[8]]

bym3.stan <- stan_model("SpatialRE/bym3.stan")
bym3.fit <- sampling(bym3.stan, data=list(N=N,y=count,E=exposure,l=l,L=mat.desgn), 
                    control = list(adapt_delta = 0.97), 
                    chains=3, warmup=7000, iter=8000, save_warmup=FALSE);
print(bym2_fit, digits=3, pars=c("beta0", "rho", "sigma", "mu[1]", "mu[2]", "mu[3]", "mu[500]", "mu[1000]", "mu[1500]", "mu[1900]", "phi[1]", "phi[2]", "phi[3]", "phi[500]", "phi[1000]", "phi[1500]", "phi[1900]", "theta[1]", "theta[2]", "theta[3]", "theta[500]", "theta[1000]", "theta[1500]", "theta[1900]"), probs=c(0.025, 0.5, 0.975));


model <- stan_model(model_code = "parameters { real<lower = 0> y; }
                           transformed parameters { real<upper = -1> z = y; }")

fit <- sampling(model)
