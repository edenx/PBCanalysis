# Stratification level construction with PBC data

split_char <- strsplit(as.character(PBC$LSOA04NM), '[[:space:]]')
lab_reg <- sapply(split_char, function(x) c(paste0(x[-length(x)], collapse=""),x[length(x)]))
lab_reg_df <- t(lab_reg) %>% as.data.frame(.) %>% 
        rename("District"=V1, "Code"=V2) %>%
        mutate("index"=1:nrow(.))

prob_grp <- lab_reg_df %>% group_by(District, add=TRUE) %>% summarise(freq=length(index))
lab_reg_df_ <- merge(lab_reg_df, prob_grp, "District") %>% mutate(prob=freq/nrow(.))
lab_reg_df_ <- lab_reg_df_[order(lab_reg_df_$index),]

sample(lab_reg_df_$index, 1, prob=lab_reg_df_$prob)