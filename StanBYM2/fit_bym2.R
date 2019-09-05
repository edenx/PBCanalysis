library(rstan)
library(rgeos)
library(rgdal)
library(spdep)
library(sf)
library(maptools)
library(Matrix)
library(INLA)

# 167650 182338 gives TRUE for the rgeos::gTouches
# 8480 8507 for the labelling in sp spatial polygon dataframe
index_mat <- rownames(PBCshp@data)
index_mat[!index_mat %in% c(8480,8507)] <- 0

plot(PBCshp)
pointLabel(coordinates(PBCshp), labels=index_mat, cex=0.7)

# go ahead and use the operation for sf object
PBC_touches <- as.list(st_touches(PBC, PBC))

# turn the list of vector into sparse adjacency matrix
n.ids <- sapply(PBC_touches, length)
vals <- unlist(PBC_touches)
vals_diff <- vals - rep(1:545, times=sapply(PBC_touches, length))
# 
n.ids_ <- rep_len(0, length(n.ids))
acc <- 0
vals_ <- c()
for(i in 1:length(n.ids)){
        acc_n <- acc + n.ids[i]

        pos_indx <- vals_diff[(acc+1):acc_n]>0
        vals_ <- c(vals_, vals[(acc+1):acc_n][pos_indx])
        
        print(vals_)
        n.ids_[i] <- sum(pos_indx)
        acc <- acc_n
}

PBC_adj <- sparseMatrix(vals_, rep(seq_along(n.ids_), n.ids_), symmetric=TRUE)
PBC_adj_Tsparse <- sparseMatrix(vals_, rep(seq_along(n.ids_), n.ids_), giveCsparse=FALSE)

# specifying the parameters for the model
count <- PBC$X
E <- exp(as.vector(scale(log(PBC$pop), center=FALSE)))
N <- nrow(PBC)
N_edges <- as.integer(sum(PBC_adj)/2)
node_1 <- as.array(PBC_adj_Tsparse@j + as.integer(1))
node_2 <- as.array(PBC_adj_Tsparse@i + as.integer(1))
design_mtx <- as.data.frame(as.matrix(PBC[c("Income", "Crime", "Environment", "Employment",
                                            "Barriers", "propmale", "Education")])[, -8])

# find the scaling factor
# The ICAR precision matrix (note! This is singular)
Q <-  Diagonal(N, rowSums(PBC_adj)) - PBC_adj
# Add a small jitter to the diagonal for numerical stability (optional but recommended)
Q_pert <- Q + Diagonal(N) * max(diag(Q)) * sqrt(.Machine$double.eps)

# inversion for sparse matrix (INLA)
Q_inv <- inla.qinv(Q_pert, constr=list(A=matrix(1,1,N), e=0))

#Compute the geometric mean of the variances, which are on the diagonal of Q.inv
scaling_factor <- exp(mean(log(diag(Q_inv))))

# get a mean estimateb for the poisson glm without ICAR
init_mod <- glm(count ~ 1 + offset(log(E)), family=poisson)
init_beta0 <- init_mod$coefficients

# Fit the stan bym2 model
options(mc.cores = parallel::detectCores())
SEED <- 727282

bym2_stanfit <- stan("bym2.stan", 
                    data=list(N=N, N_edges=N_edges, node1=node_1, node2=node_2, y=count, 
                              # x=design_mtx, 
                              E=E, scaling_factor=scaling_factor, mu_0=init_beta0), 
                    control=list(max_treedepth=20, adapt_delta=0.97), thin=10,
                    chains=4, warmup=6000, iter=12000, save_warmup=FALSE);

shinystan::launch_shinystan(bym2_stanfit)

print(bym2_stanfit, pars=c("lp__", "beta0", "rho", "sigma", "log_precision", 
                           "logit_rho", "mu[5]", "phi[5]", "theta[5]"), 
      probs=c(0.025, 0.5, 0.975))


