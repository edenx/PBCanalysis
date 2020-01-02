library(tidyverse)
library(sf)
library(INLA)

design_mtx <- as.data.frame(as.matrix(PBC[c("Income", "Crime", "Environment", "Employment",
                                            "Barriers", "propmale", "Education")])[, -8])

M <- spdep::poly2nb(PBC) %>% 
        spdep::nb2mat(style = "B", zero.policy = TRUE)

colnames(M) <- rownames(M)

#' Scaled precision matrix for BYM2


f_inla <- inla(PBC$X ~ f(seq_along(M[1,]), model = "bym2", graph = M),
               family = "binomial", data = PBC, Ntrials = nrow(PBC),
               control.predictor=list(link=NA),
               control.compute=list(config = TRUE))

#' Sample from posterior distribution
smp_inla <- inla.posterior.sample(1000, f_inla, intern = TRUE)

#' Simulate logit prevalence and output prevalence
mu_rho_inla <- lapply(smp_inla, "[[", "latent") %>%
        vapply("[", numeric(nrow(PBC)), seq_len(nrow(PBC)))

plhiv_inla <- plogis(mu_rho_inla) * PBC$pop
rho_out_inla <- (matrix(1, 1, nrow(M)) %*% plhiv_inla) / as.numeric(matrix(1, 1, nrow(M)) %*% PBC$pop)

sim_inla <- list(mu_rho = mu_rho_inla,
                 rho_out = rho_out_inla)


f_inla$summary.linear.predictor %>%
  transmute(mean,
            sd,
            median = `0.5quant`,
            lower = `0.025quant`,
            upper = `0.975quant`)
