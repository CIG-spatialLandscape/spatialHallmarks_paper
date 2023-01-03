library(dplyr)
df <- readxl::read_xlsx("Downloads/Rstudiotumor(1).xlsx")
df$Others[20] <- "100"
colnames(df)[1] <- "clusters"
head(df)


df_filtered <- df %>% filter(!is.na(`Tumor %`)) %>% select(clusters,`Tumor %`, `Stroma %`, `Necrosis %`)
library(reshape2)
#df_filtered <- df %>% filter(!is.na(`Tumor %`)) %>% select(clusters,`Tumor %`, `Stroma %`, `Necrosis %`) %>% melt()

library(ggplot2)
library(ggpubr)
p1 <- ggplot(df_filtered, aes(x=clusters, y=`Tumor %`, fill=clusters)) + geom_boxplot() +
  scale_fill_manual(values = c("cluster 5" = "lightgoldenrod1", "cluster 4" = "lightgoldenrod3", 
                               "cluster 3" = "lightpink2", "cluster 2"= "orchid3", "cluster 1"= "orchid4")) + 
  theme_pubr(base_size = 20) + theme(legend.position = "none") + labs(x="", y="Cancer cells %") + scale_x_discrete(labels=1:5)

p2 <- ggplot(df_filtered, aes(x=clusters, y=`Stroma %`, fill=clusters)) + geom_boxplot() +
  scale_fill_manual(values = c("cluster 5" = "lightgoldenrod1", "cluster 4" = "lightgoldenrod3", 
                               "cluster 3" = "lightpink2", "cluster 2"= "orchid3", "cluster 1"= "orchid4")) + 
  theme_pubr(base_size = 20) + theme(legend.position = "none") + labs(x="", y="Stroma cells %") + scale_x_discrete(labels=1:5)

p1/p2


df$`Normal %` <- readr::parse_number(df$Others)
df$`Normal %`[is.na(df$`Normal %`)] <- 0
df_filtered <- df %>% filter(!is.na(`Tumor %`)) %>% select(clusters,`Tumor %`, `Stroma %`, `Necrosis %`, `Normal %`) %>% melt()
df_filtered$variable <- as.character(df_filtered$variable)
df_filtered$variable[df_filtered$variable == "Normal %"] <- "Normal cells"
df_filtered$variable[df_filtered$variable == "Tumor %"] <- "Tumor cells"
ggpaired(df_filtered[!df_filtered$variable %in% c("Stroma %", "Necrosis %"),], x = "variable", y = "value",
         color = "variable", palette = "jco", 
         line.color = "gray", line.size = 0.4,
         facet.by = "clusters", short.panel.labs = T, ggtheme = theme_pubr())  + labs(x="", y="%", col="Type") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + facet_wrap(~clusters, nrow = 1)
######################################3

ggplot(df_filtered, aes(x=clusters, y=`Necrosis %`, fill=clusters)) + geom_boxplot() +
  scale_fill_manual(values = c("cluster 5" = "lightgoldenrod1", "cluster 4" = "lightgoldenrod3", 
                               "cluster 3" = "lightpink2", "cluster 2"= "orchid3", "cluster 1"= "orchid4")) + 
  theme_classic()
ggplot(df_filtered[df_filtered$variable != "Necrosis %",], aes(x=clusters, y=value, fill=variable)) + geom_boxplot() + geom_jitter()
  theme_classic() 

df_filtered$new <- df_filtered$`Tumor %`/df_filtered$`Stroma %`
df_filtered$new[df_filtered$new == Inf] <- 20
df_filtered$new[df_filtered$new == NaN] <- 0
ggplot(df_filtered, aes(x=clusters, y=new, fill=clusters)) + geom_boxplot() +
  scale_fill_manual(values = c("cluster 5" = "lightgoldenrod1", "cluster 4" = "lightgoldenrod3", 
                               "cluster 3" = "lightpink2", "cluster 2"= "orchid3", "cluster 1"= "orchid4")) + 
  theme_classic() + geom_hline(yintercept = 1)



library(ggpubr)
p <- ggpaired(df_filtered[df_filtered$variable != "Stroma %",], x = "variable", y = "value",
              color = "variable", palette = "jco", 
              line.color = "gray", line.size = 0.4,
              facet.by = "clusters", short.panel.labs = FALSE)



###

v <- readRDS("Desktop/server/ClusterPlot/data/v.rds")
v <- v$Kidney3

hires<- readRDS("Desktop/server/ClusterPlot/data/hires/Kidney3.rds")


tmp <- v$imagecol
v$imagecol <- v$imagerow
v$imagerow <- tmp


tmp <- ref_v$imagecol
ref_v$imagecol <- ref_v$imagerow
ref_v$imagerow <- tmp
