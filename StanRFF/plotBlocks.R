
lis_blocks <- list()
sf_geom <- PBC[, c("geometry", "X")] %>% 
        rename(count=X)
filt_index <- c(1,0,1,0,1,0,1,1,1,1)==1

for(i in seq(length(fin_index)-1)[filt_index]){
        
        # png(paste0("Block", fin_index[i], "_", fin_index[i+1], ".png"), width=500)
        
        sf_pred <- PBC[seq(fin_index[i], fin_index[i+1]), c("geometry", "X")] %>% 
                # mutate(pred=pred, normpop=lis_wcentr) %>% 
                rename(count=X)
        
        pop_ <- ggplot() + 
                geom_sf(data=sf_geom["count"], lwd=0.05) +
                geom_sf(data=sf_pred["count"], aes(fill=count), lwd=0.05) + 
                scale_fill_distiller(palette="BuPu", direction=1) +
                theme_linedraw() +
                labs(x="Northing", y="Easting")
        print(pop_)
        # dev.off()
}
