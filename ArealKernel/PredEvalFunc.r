# --------------------------------- Evaluation Metrics ------------------------------------------
# Bias
bias <- function(pred) mean(colMeans(pred-count))

# RMSE
rmse <- function(pred) sqrt(mean(colMeans((pred-count)^2)))

# Credible Interval 95%
cp95 <- function(pred){
        ci95 <- apply(pred, 2, quantile, c(0.975, 0.025))
        mean(ci95[1,]>=count & ci95[2, ]<=count)
}

# prediction: histogram
plot_hist <- function(pred, alpha){
        col1 <- adjustcolor(2, alpha.f = 0.3)
        col2 <- adjustcolor(4, alpha.f = 0.3)
        hist(pred, col=col1, 
             main=paste0("Data vs. Predicted Mean Count with ls=", alpha))
        hist(PBC_df$X, col=col2, add=TRUE)
        legend("topright", legend=c("Predicted Mean", "Count"),
               fill=c(col1, col2), cex=0.8)
}