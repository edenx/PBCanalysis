# ----------------------------- Producing aggregated Kernel -------------------------------------

# for every intersected cells with polygon, normalise the population density per polygon
check.0 <- function(x) {
        w_val <- x[,"value"]*x[,"weight"]
        if(any(x[,"value"] != 0)){
                cbind(x, W=x[,"value"]/sum(w_val))}
        else{cbind(x, W=x[,"value"])}
}

# distance matrix between cells within two regions A, B
vectorised_pdist <- function(A, B){
        an <- apply(A, 1, function(rvec) crossprod(rvec,rvec))
        bn <- apply(B, 1, function(rvec) crossprod(rvec,rvec))
        
        m <- nrow(A)
        n <- nrow(B)
        
        tmp <- matrix(rep(an, n), nrow=m) 
        tmp <- tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
        sqrt(tmp - 2*tcrossprod(A,B))
}

# alpha lengthscale (maybe a bit confusing, change afterwards)
# nu any in c(5, 3, 1, "inf")
Matern.ker <- function(coord1, coord2, alpha, nu){
        # dim: nrow(coord1) x nrow(coord2)
        d <- vectorised_pdist(coord1, coord2)
        
        if(!nu %in% c(1, 3, 5, "inf")){
                errorCondition("choose appropriate nu in k52, k32, k12, inf")
        }else if(nu==5){
                (1 + sqrt(5)*d/alpha + 5*d^2/(3*alpha^2)) * exp(-sqrt(5)*d/alpha)
        }else if(nu==3){
                (1 + sqrt(3)*d/alpha) * exp(-sqrt(3)*d/alpha)
        }else if(nu==1){
                exp(-d/alpha)
        }else{
                exp(-d^2/alpha^2)
        }
}

Region.ker <- function(wcentr1, wcentr2, alpha, nu=5){
        centr1 <- as.matrix(wcentr1[,c("x", "y")])
        centr2 <- as.matrix(wcentr2[,c("x", "y")])
        
        # population density empirical distribution
        weight1 <-  as.matrix(wcentr1[,c("W")]) * as.matrix(wcentr1[,c("weight")])
        weight2 <-  as.matrix(wcentr2[,c("W")]) * as.matrix(wcentr2[,c("weight")])
        
        base.ker <- Matern.ker(centr1, centr2, alpha, nu)
        
        # weighted average over region with normalised population density
        t(weight1) %*% base.ker %*% weight2
}

# Full kernel
sim.regker <- function(alpha, lis_region, nu=5, plot=FALSE){
        # parallel the outer for loop
        require(parallel)
        no_cores <- detectCores() - 1
        cl <- makeCluster(no_cores, type="FORK")
        
        Kernel <- parSapply(cl, lis_region, function(region) 
                sapply(lis_region, Region.ker, region, alpha, nu))
        stopCluster(cl)
        
        if(plot){
                print(lattice::levelplot(Kernel))
        }
        
        return(Kernel)
}

rff.region <- function(wcentr, Omega, m, alpha){
        centr <- as.matrix(wcentr[,c("x", "y")])
        # population density empirical distribution
        w <- as.matrix(wcentr[,c("W")])
        h <- as.matrix(wcentr[,c("weight")])
        
        # Projection - combine data with sample frequencies
        proj <- centr %*% t(Omega) 
        
        # Fourier feature for a given area
        # print(w*h)
        phi <- sqrt(1/m) * colSums(t(w*h) %*% cbind(cos(proj/alpha), sin(proj/alpha)))
        
        return(phi)
}

# RFF approximated kernel (Cholesky decomposition)
# feature sampling: using the inversion -- quasi monte carlo integration;
# Matern 5-2 kernel -- student 5-2 sdf
sim.rff <- function(alpha, lis_region, m=100, nu=5, plot=FALSE){
        Omega <- qt(halton(m,2), nu)
        Phi <- t(sapply(lis_region, rff.region, Omega, m, alpha))
        
        Kernel <- Phi %*% t(Phi)
        if(plot){
                print(lattice::levelplot(Kernel))
                return(as.matrix(Kernel))
        }
        
        return(as.matrix(Phi))
}