# following your example code, I used exponential covariance function, i.e matern 0.5 and it worked
u = seq(0,1,.01)
exact.cov = exp(-u/0.25) # using length scale of 0.25
plot(u,exact.cov,ty="l")

# now let us build an approximation:
freq = rt(100,1.0)/.25  # 1 degree of freedom work for matern 0.5
Kx = sqrt(1/(length(freq))) * cbind(cos(matrix(u)  %*% freq),
                                    sin(matrix(u)  %*% freq))
Kapprox = Kx %*% t(Kx)
lines(u,Kapprox[1,],col="blue")

# and a better approximation
freq = rt(10000,1.0)/.25
Kx = sqrt(1/(length(freq))) * cbind(cos(matrix(u)  %*% freq),
                                    sin(matrix(u)  %*% freq))
Kapprox = Kx %*% t(Kx)
lines(u,Kapprox[1,],col="red")

################################## Now experimenting in a geostatistical setting #######################
# I want to see what the covariance matrix will look like ###########
n <- 100   #number of point
# generate the coordinates ######
xy <- as.matrix(expand.grid(x=seq(0, 1, length.out = sqrt(n)), y=seq(0, 1, length.out = sqrt(n)))) # create the coordinates
u <- as.matrix(dist(xy)) # create the distance matrix

k52 = function(d,lengthscale) {
  (1 + sqrt(5)*d/lengthscale + 5*d^2/(3*lengthscale^2)) * exp(-sqrt(5)*d/lengthscale)
}
k32 = function(d,lengthscale) {
  (1 + sqrt(3)*d/lengthscale) * exp(-sqrt(3)*d/lengthscale)
}
k12 = function(d,lengthscale) {
  exp(-d/lengthscale)
}

u1 <- as.matrix(dist(xy[,1])) # create the distance matrix
u2 <- as.matrix(dist(xy[,2])) # create the distance matrix
exact.cov = k52(u,1/600)
lattice::levelplot(exact.cov)

#exact.cov = k12(u1,1)*k12(u2,1)
#exact.cov = k52(u,1)
#exact.cov <- exp(-.5*u^2) # exp(-u)
# 
# build an approximation similar to the 1 dimension
m <- 1000 # number of frequencies
freq <- stats::rt(m*2, df=5)   # using length scale of 0.25
# create a matrix of the frequencies with row corresponding to the freq_1 \ldots freq_m
mat.freq <- matrix(freq, nrow=m, ncol=2, byrow=FALSE)
PHI <- cbind(cos(xy%*%t(mat.freq)), sin(xy%*%t(mat.freq)))
CovAprox <- (PHI%*%t(PHI))/m
lattice::levelplot(CovAprox)

library(mvnfast)
# Building a better approximation by increasing the number of frequency
m2 <- 100000 # number of frequencies
#freq2 <- rnorm(m2*2) #stats::rt(m2*2, df=5)    # using length scale of 0.25
freq2 <- rmvt(n=m2, c(0,0), diag(1,2), df=5) # stats::rt(m2*2, df=1)    # using length scale of 0.25
# create a matrix of the frequencies with row corresponding to the freq_1 \ldots freq_m
mat.freq2 <- matrix(freq2, nrow=m2) #, ncol=2, byrow=FALSE)
PHI2 <- cbind(cos(scale(xy)%*%t(mat.freq2)), sin(scale(xy)%*%t(mat.freq2))) 
CovAprox2 <- (PHI2%*%t(PHI2))/m2

lattice::levelplot(CovAprox2)
# 
# #################
# plot(u[1,],exact.cov[1,],ty="p")
# points(u[1,],CovAprox[1,],col="blue")
# points(u[1,],CovAprox2[1,],col="red")
#######
par(mfrow=c(2,2))
ii = 1:sqrt(n)
plot(u[1,ii],exact.cov[1,ii],ty="l")
lines(u[1,ii],CovAprox[1,ii],col="blue")
lines(u[1,ii],CovAprox2[1,ii],col="red")

ii = ii + sqrt(n)
plot(u[1,ii],exact.cov[1,ii],ty="l")
lines(u[1,ii],CovAprox[1,ii],col="blue")
lines(u[1,ii],CovAprox2[1,ii],col="red")


ii = ii + 2*sqrt(n)
plot(u[1,ii],exact.cov[1,ii],ty="l")
lines(u[1,ii],CovAprox[1,ii],col="blue")
lines(u[1,ii],CovAprox2[1,ii],col="red")

ii = ii + 2*sqrt(n)
plot(u[1,ii],exact.cov[1,ii],ty="l")
lines(u[1,ii],CovAprox[1,ii],col="blue")
lines(u[1,ii],CovAprox2[1,ii],col="red")


