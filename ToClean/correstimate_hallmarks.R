df_all

c <- c()
df_cor <- data.frame()
for (sample in unique(df_all$sample)) {
  df_cor <- rbind(df_cor, c(cor(df_all[df_all$sample==sample, "H1"], df_all[df_all$sample==sample, "estimate"]),
                  cor(df_all[df_all$sample==sample, "H2"], df_all[df_all$sample==sample, "estimate"]),
                  cor(df_all[df_all$sample==sample, "H3"], df_all[df_all$sample==sample, "estimate"]),
                  cor(df_all[df_all$sample==sample, "H4"], df_all[df_all$sample==sample, "estimate"]),
                  cor(df_all[df_all$sample==sample, "H5"], df_all[df_all$sample==sample, "estimate"]),
                  cor(df_all[df_all$sample==sample, "H6"], df_all[df_all$sample==sample, "estimate"]),
                  cor(df_all[df_all$sample==sample, "H7"], df_all[df_all$sample==sample, "estimate"]),
                  cor(df_all[df_all$sample==sample, "H8"], df_all[df_all$sample==sample, "estimate"]),
                  cor(df_all[df_all$sample==sample, "H9"], df_all[df_all$sample==sample, "estimate"]),
                  cor(df_all[df_all$sample==sample, "H10"], df_all[df_all$sample==sample, "estimate"]),
                  cor(df_all[df_all$sample==sample, "H11"], df_all[df_all$sample==sample, "estimate"]),
                  cor(df_all[df_all$sample==sample, "H12"], df_all[df_all$sample==sample, "estimate"]),
                  cor(df_all[df_all$sample==sample, "H13"], df_all[df_all$sample==sample, "estimate"]))
                  )
}
colnames(df_cor) <- paste0("H", 1:13)
library(reshape2)
long_df_cor <- melt(df_cor)
long_df_cor$specific <- "TME"
long_df_cor$specific[long_df_cor$variable %in% c("H2", "H4", "H8", "H10", "H11", "H12")] <- "Cancer"
rownames(palette) <- NULL
ggplot(long_df_cor, aes(x=variable, y=value, fill=variable)) + geom_boxplot() + 
  facet_grid(~specific, scales = "free") + theme_classic() + labs(x="Hallmarks", y="Correlation (Hallmark activity & ESTIMATE)") + 
  scale_fill_manual(values = palette[,1]) + coord_flip()
