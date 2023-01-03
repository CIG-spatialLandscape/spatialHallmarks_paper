### Circular barplot ###

# library
library(tidyverse)
library(dplyr)
library(Seurat)

# Define CancerType groups
lm_mat.basal.sub$CancerType <- NA
lm_mat.basal.sub$CancerType[lm_mat.basal.sub$Sample %in% c("Acinar", "IC")] <- "Prostate"
lm_mat.basal.sub$CancerType[lm_mat.basal.sub$Sample %in% c("Breast", "BreastA", "Ductal", "DuctalFFPE", "TNBCA")] <- "Breast"
lm_mat.basal.sub$CancerType[lm_mat.basal.sub$Sample %in% c("C20", "C21", "C34", "C51", "C7")] <- "Kidney"
lm_mat.basal.sub$CancerType[lm_mat.basal.sub$Sample %in% c("cHC1T", "HCC1T", "HCC2T", "HCC5D", "ICC1L")] <- "Liver"
lm_mat.basal.sub$CancerType[lm_mat.basal.sub$Sample %in% c("Colorectal", "CRC1", "CRC2", "Intestine")] <- "Colorectal"
lm_mat.basal.sub$CancerType[lm_mat.basal.sub$Sample %in% c("DU2", "DU3", "DU12", "DU8", "DU13")] <- "Bladder"
lm_mat.basal.sub$CancerType[lm_mat.basal.sub$Sample %in% c("OV4A", "OVD1", "Ovarian", "CUP295", "OVFFPE")] <- "Ovarian"
lm_mat.basal.sub$CancerType[lm_mat.basal.sub$Sample %in% c("P259_H2A2", "P264", "P270", "P288", "P306")] <- "PDAC"
lm_mat.basal.sub$CancerType[lm_mat.basal.sub$Sample %in% c("Glioblastoma", "UKF242T", "UKF260T", "UKF269T", "UKF275T")] <- "Glioblastoma"



### Major dendogram ###

library(dendextend)
library(circlize)

df <- data.frame(matrix(nrow = 13, ncol = 41))
colnames(df) <- unique(lm_mat.basal.sub$Sample)
rownames(df) <- paste0("H",1:13)
for (h in paste0("H",1:13)){
  for (sample in unique(lm_mat.basal.sub$Sample)) {
    df[h, sample] <- lm_mat.basal.sub[lm_mat.basal.sub$Sample==sample & lm_mat.basal.sub$Hallmark == h, "diff"]
  }
}

# Distance matrix
d <- dist(df)
# Hierarchical clustering dendrogram

hc <- hclust(d, method = "average")
hallmark_order <- hc$labels[hc$order]

hc <- as.dendrogram(hc)

# Circular dendrogram

png("Desktop/barplot/dendrograms/hallmarks_dendrogram.png")
circlize_dendrogram(hc %>% set("branches_lwd", 7),
                    labels_track_height = NA,
                    dend_track_height = 0.5, width=50)
dev.off()



### Individual dendrogram ###
tmp <- lm_mat.basal.sub[lm_mat.basal.sub$Hallmark == "H12", ] %>% group_by(CancerType) %>% summarise(mean = mean(diff), sd = sd(diff))
rownames(tmp) <- tmp$CancerType
# Distance matrix
d <- dist(tmp)
# Hierarchical clustering dendrogram
hc <- hclust(d, method = "average")
tmp$CancerType[hc$order]

lm_mat.basal.sub <- lm_mat.basal.sub[order(lm_mat.basal.sub$CancerType),]

lm_mat.basal.sub 

lm_mat.basal.sub %>%
  arrange(factor(CancerType, levels = tmp$CancerType[hc$order]))



new_df <- c()
for (hallmark in paste0("H", 1:13)) {
  tmp <- lm_mat.basal.sub[lm_mat.basal.sub$Hallmark == hallmark, ] %>% group_by(CancerType) %>% summarise(mean = mean(diff), sd = sd(diff))
  rownames(tmp) <- tmp$CancerType
  # Distance matrix
  d <- dist(tmp)
  # Hierarchical clustering dendrogram
  hc <- hclust(d, method = "average")
  png(paste0("Desktop/barplot/dendrograms/", hallmark, "_dendrogram.png"))
  plot(hc,hang = -1, lwd = 7)
  dev.off()
  new_df <- rbind(new_df, lm_mat.basal.sub[lm_mat.basal.sub$Hallmark == hallmark, ] %>%
                    arrange(factor(CancerType, levels = tmp$CancerType[hc$order])))  
  
}


data <- new_df
data$Hallmark <- factor(data$Hallmark, levels = hallmark_order)

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 5
to_add <- data.frame( matrix(NA, empty_bar*13, ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$Hallmark <- rep(paste0("H", 1:13), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(Hallmark)
data$id <- seq(1, nrow(data))

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(Hallmark) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))
tmp <- data[!is.na(data),] %>% group_by(Hallmark) %>% summarise(mean = mean(diff), sd = sd(diff))
base_data$mean <- tmp$mean[1:13]

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

data$model <- factor(data$model, levels = c("increase", "basal"))
# Make the plot
p <- ggplot(data, aes(x=as.factor(id), y=diff, fill=model)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", width = 1) +
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0.5, xend = start, yend = 0.5), colour = "grey", alpha=0.8, size=0.2 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.2, xend = start, yend = 0.2), colour = "grey", alpha=0.8, size=0.2 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.1, xend = start, yend = 0.1), colour = "grey", alpha=0.8, size=0.2 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),3), y = c(0.1, 0.2, 0.5), label = c("0.1", "0.2", "0.5") , color="grey", size=3 , angle=0, fontface="bold", hjust=0.9) +
  scale_fill_manual(values = DiscretePalette(9)) + 
  ylim(-1.2,1) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-3,4), "cm")
  ) + 
  #annotate("text", y = -0.5, x = 0, label = "R-squared gain by Neighborhood", size=4.5) + 
  coord_polar() + 
  geom_segment(data=base_data, aes(x = start, y = mean, xend = end, yend = mean), colour = "black", alpha=1, size=0.4 , inherit.aes = FALSE )   + 
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -0.01, xend = end, yend = -0.01), colour = "black", alpha=0.8, size=0.5 , inherit.aes = FALSE )  
#geom_text(data=base_data, aes(x = title, y = -0.05, label=Hallmark), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

p


ggsave("Desktop/barplot/barplot.png", bg = "white", width = 9, height = 9)  



###
data$model <- "increase"
a <- data
a$diff <- a$Rsquared
a$model <- "basal"
data <- rbind(data, a)

tmp <- data[data$model == "basal", ]
tmp$total <- tmp$diff + tmp$Rsquared
