require("MASS")
require("dplyr")
require("ggplot2")
require("SDALGCP")

set.seed(727261)

# Computational Grid of the PBC data
Grid_PBC <- as.data.frame(coordinates(pop_den))
n <- nrow(Grid_PBC)

# Simulate covariate variables
x1 <- rnorm(n) %>% round(2)
x2 <- rnorm(n) %>% round(2)

# Simulate spatial field (exp kernel)
distance <- as.matrix(dist(Grid_PBC))
Sigma <- 0.4*exp(-0.1*distance)
Cholesk_Sigma <- chol(Sigma)
#       Grid too large to process
omega <- t(Cholesk_Sigma) %*% MASS::mvrnorm(mu=rep(0,n), Sigma=diag(n))
# omega <- MASS::mvrnorm(mu=rep(0,n), Sigma=0.4*exp(-0.1*distance))

# Poisson distributed given random intensity eta
eta <- x1 + x2 + omega
dat <- Grid_PBC %>% mutate(y=rpois(n,exp(eta)))

# Spatial surface
ggplot(dat %>% mutate(omega=omega),
       aes(x=x.easting,y=x.northing)) +
        geom_raster(aes(fill=omega)) +
        theme_bw()

# Normal Data
ggplot(dat, aes(x=x.easting,y=x.northing,fill=y)) +
        geom_raster() +
        theme_bw()


#       Cholesky decomposition robust for this Covariance matrix
table((t(Cholesk_Sigma) %*% Cholesk_Sigma - 0.4*exp(-0.1*distance)) < 10^(-6))
