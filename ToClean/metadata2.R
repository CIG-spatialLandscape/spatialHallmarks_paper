files <- list.files("Desktop/IJC/datasets/IGTP/figuresPaper/SCD/")
strsplit(files, split = ".", fixed = T)[[1]][1]


tmp <- sapply(files, function(file) {
  read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/SCD/", file))
})
df <- do.call(rbind, tmp)
rownames(df) <- sapply(files, function(files){strsplit(files, split = ".", fixed = T)[[1]][1]})
colnames(df) <- "SCD"

df_tmp <- df_all %>% select(c(sample, clusters))
df_tmp$compartment <- "Cancer"
df_tmp$compartment[df_tmp$clusters %in% c(4,5)] <- "TME"
df_tmp$compartment[df_tmp$clusters %in% c(3)] <- "Buffer"

df <- cbind(df_tmp %>% group_by(sample, compartment) %>% count() %>% 
  group_by(sample) %>% mutate(perc=n/sum(n)) %>% select(c(compartment, perc)) %>%
  spread(compartment, perc) %>% ungroup() %>% select(Buffer, Cancer, TME) %>% as.data.frame() , df)
write.table(df, "Desktop/df_meta.txt", sep = "\t", quote = F)
