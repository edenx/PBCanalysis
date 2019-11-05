Constr_Cmat <- function(Phi){
        list(Phi %*% t(Phi), Diagonal(nrow(Phi), 1))
}
Cmat_ <- lapply(lis_Phi.null, Constr_Cmat)

PBC_df.desn <- PBC_df %>% 
        mutate(ID=1:nrow(.), count=.$X) %>%
        dplyr::select(ID, count, Income, Crime, Environment, Employment, 
                      Barriers, propmale, Education)
# offset
pop <- scale(PBC_df$pop, center=FALSE)
# model formula
lmod <- lm(count ~ .-ID, PBC_df.desn)
form <- formula(lmod)

prec.prior <- list(prec = list(param = c(0.001, 0.001)))

inla_reg.meff_ <- function(Phi, mat_desn, form, Cmat, offset_){
        fit <- inla(
                update(form,. ~. + f(ID, model="generic3", Cmatrix=Cmat_)),
                data=mat_desn, family = "poisson", E=offset_,
                control.predictor = list(compute=TRUE))
        
        fv <- fit$summary.fitted.values
        mlik <- fit$mlik[1]
        
        out <- list(fv=fv, mlik=mlik, fit=fit)
        
        return(out)
}