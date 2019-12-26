ibrary(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

N = 1000
x = runif(N,0,2)
x = x[order(x)]
y = sin(2*pi*x) + rnorm(N,0,.3)

k=10
data = list(y=y,x=x,n=length(y),k=k,bw=8,omega=rnorm(k))

# fixed bandwidth
m1 = stan_model("rff1.stan")
fit = sampling(m1, data=data, iter=100, warmup=50, chains=4)
out = extract(fit) 

png("figure.png")
par(mfrow=c(1,2))
plot(x,y,pch=18,col="gray",main="Fixed lengthscale")
lines(x,colMeans(out$fhat),col="blue")

# learn the bandwidth
m2 = stan_model("rff2.stan")
fit = sampling(m2, data=data, iter=200, warmup=100, chains=4)
out = extract(fit) 
print(fit,c("bw","lp__"))

plot(x,y,pch=18,col="gray",main="Learnt lengthscale")
lines(x,colMeans(out$fhat),col="red")
dev.off()
