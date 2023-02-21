library(dplyr)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(stringr)
library(tibble)

source("Desktop/IJC/git/spatialHallmarks_paper/utils/SamplesMetadata.R")
df <- read.table("Desktop/toClean/Hallmark_autocorrelation.txt", sep = "\t")

colnames(df)[2:14] <- hallmark_names
#ridge plot for each hallmark
ggplot(df %>% reshape2::melt(id.vars="sample"), aes(y=variable, x=value, fill=variable)) + geom_density_ridges() + xlim(c(0,1)) +
  scale_fill_manual(values = do.call(c, color_codes))+ theme_ridges(center_axis_labels = T) + theme(legend.position = "none") + labs(x="Moran's I", y="")

#Moran's I within each compartment
samples <- df$sample
df_all_samples <- c()
for (sample_id in samples) {
  df_whole <- df %>% filter(sample==sample_id) %>% remove_rownames %>% column_to_rownames(var="sample") %>% t()  %>% as.data.frame() 
  colnames(df_whole) <- c("cor")
  df_whole$tissue <- "Whole"
  df_whole$signature <- rownames(df_whole)
  
  
  #Neoplastic 
  df_neoplastic <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/moran/Neoplastic/", sample_id,"_HallmarkMoran.txt"), sep = "\t")
  df_neoplastic <- df_neoplastic %>% filter(!hallmark %in% paste0("H", c(1,3,5,6,7,9,13)))

  df_neoplastic <- df_neoplastic %>% select(-hallmark)
  df_neoplastic$tissue <- "Neoplastic"
  df_neoplastic$signature <- rownames(df_neoplastic)
  
  #TME 
  
  df_TME <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/moran/TME/", sample_id,"_HallmarkMoran.txt"), sep = "\t")
  df_TME <- df_TME %>% filter(hallmark %in% paste0("H", c(1,3,5,6,7,9,13)))

  df_TME <- df_TME %>% select(-hallmark)
  df_TME$tissue <- "TME"
  df_TME$signature <- rownames(df_TME)
  
  df_all <- rbind(df_whole, df_neoplastic, df_TME)
  df_all$id <- paste0(df_all$tissue, "_", df_all$signature2)
  df_all_samples <- rbind(df_all_samples, df_all)
}
#boxplot 
df_all_samples$tissue <- factor(df_all_samples$tissue, levels = c("Whole","Neoplastic", "TME"))
#df_all_samples %>% group_by(tissue) %>% summarise(m=mean(cor))
ggplot(df_all_samples, aes(x=tissue,y=cor, fill=tissue))  + 
  geom_violin() + geom_boxplot(width=0.1, outlier.size = 0.4, alpha=0.2, ) + ylim(c(0,1.05)) + 
  stat_compare_means(size=6, family="Times New Roman") + theme_classic() +
  labs(y="Moran's I", x="Tissue") + 
  ggtitle("Spatial autocorrelation of Hallmark activities") +
  scale_fill_manual(values = c("gray","orchid4", "lightgoldenrod1")) + 
  theme(axis.text = element_text(size=18), #axis.text.x = element_text( angle = -45, hjust = 0),
        axis.title = element_text(size=14, face="bold"),
        legend.position = "none", plot.title = element_text(size=18,hjust=0.5,face="bold"))
ggsave("Desktop/IJC/datasets/IGTP/figuresPaper/Final/moran_all.png", width = 7, height = 5)

