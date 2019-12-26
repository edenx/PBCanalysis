require("MASS")
require("dplyr")
require("ggplot2")
require("SDALGCP")
set.seed(727261)

### Prepare the input of the model
data <- as.data.frame(PBCshp@data)  #get the data

### Write the formula of the model
FORM <- X ~ propmale + Income + Employment + Education + Barriers + Crime +
        Environment +  offset(log(pop))

# FORM <- X ~ 1

### set the discretised phi
phi <- seq(500, 1700, length.out = 20)

#### get the initial parameter
model <- glm(formula=FORM, family="poisson", data=data)
beta.start <-coef(model)
sigma2.start <- mean(model$residuals^2)
phi.start <- median(phi)
par0 <- c(beta.start, sigma2.start, phi.start)

# setup the control arguments for the MCMC
n <- nrow(data)
h <- 1.65/(n^(1/6))
control.mcmc <- controlmcmcSDA(n.sim = 10000, burnin = 2000,
                               thin= 8, h=h, c1.h = 0.01, c2.h = 1e-04)
###Run the model
# pop_den <- values(SDALGCP::pop_den)
data(pop_den)
pop_den[is.na(pop_den[])] <- 0
my_est <- SDALGCPMCML(formula=FORM, data=data, my_shp=PBCshp, delta=100, 
                      phi=phi, method=1, 
                      pop_shp=pop_den,
                      weighted=TRUE, 
                      plot=TRUE, par0=NULL, control.mcmc=control.mcmc)
summary(my_est)
Con_pred_2 <- SDALGCPPred(para_est=my_est,  cellsize=300, continuous=TRUE)

#to plot the spatially continuous relative risk
plot(Con_pred_2, type="relrisk")
#to plot the incidence
plot(Con_pred_2, type="incidence", continuous=FALSE)
#to plot the exceedance probability of the relative risk
plot(Con_pred, type="relrisk", thresholds= 2)
#to plot the exceedance probability of the incidence
plot(Con_pred, type="incidence", continuous=FALSE, thresholds= 0.001)