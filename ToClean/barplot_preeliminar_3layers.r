### Circular barplot ###

# library
library(tidyverse)
library(dplyr)
library(Seurat)




#merge in one table all important information
lm_mat.basal$diff <- lm_mat.1$Rsquared - lm_mat.basal$Rsquared
lm_mat.basal$RsquaredCancer <- lm_mat.basal.Cancer$Rsquared
lm_mat.basal$diffCancer <- lm_mat.1.Cancer$Rsquared - lm_mat.basal.Cancer$Rsquared
lm_mat.basal$RsquaredTME <- lm_mat.basal.TME$Rsquared
lm_mat.basal$diffTME <- lm_mat.1.TME$Rsquared - lm_mat.basal.TME$Rsquared
lm_mat.basal$RsquaredBuffer <- lm_mat.basal.Buffer$Rsquared
lm_mat.basal$diffBuffer <- lm_mat.1.Buffer$Rsquared - lm_mat.basal.Buffer$Rsquared

########### RUN from here !!!!!!!!! #########################3

lm_mat.basal <- read.table("")



# Define CancerType groups

lm_mat.basal$CancerType <- NA
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("Acinar", "IC", "PC1", "PC2")] <- "Prostate"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("Breast", "BreastA", "Ductal", "DuctalFFPE", "TNBCA", "M1", "M2", "M3", "M4")] <- "Breast"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("C20", "C21", "C34", "C51", "C7")] <- "Kidney"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("cHC1T", "HCC1T", "HCC2T", "HCC5D", "ICC1L")] <- "Liver"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("Colorectal", "CRC1", "CRC2", "Intestine", "Co1", "Co2", "Co3", "Co4")] <- "Colorectal"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("DU2", "DU3", "DU12", "DU8", "DU13")] <- "Bladder"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("OV4A", "OVD1", "Ovarian", "CUP295", "OVFFPE")] <- "Ovarian"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("P259_H2A2", "P264", "P270", "P288", "P306")] <- "Pancreas"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("Glioblastoma", "UKF242T", "UKF260T", "UKF269T", "UKF275T")] <- "Glioblastoma"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("P1","P3", "P4", "P5", "P6", "P8")] <- "Lung"

### Major dendogram ###

#sort hallmarks
library(dendextend)
library(circlize)

df <- data.frame(matrix(nrow = 13, ncol = 41))
colnames(df) <- unique(lm_mat.basal$Sample)
rownames(df) <- paste0("H",1:13)
for (h in paste0("H",1:13)) {
  for (sample in unique(lm_mat.basal$Sample)) {
    df[h, sample] <- lm_mat.basal[lm_mat.basal$Sample==sample & lm_mat.basal$Hallmark == h, "diff"]
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

data <- lm_mat.basal[order(lm_mat.basal$CancerType),]
data$Hallmark <- factor(data$Hallmark, levels = hallmark_order)

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 10
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


# get locations for each cancer group to put color lines
GroupCol <- c()
for (h in paste0("H", 1:13)) {
  base_data2 <- data[data$Hallmark==h & !is.na(data$Sample),] %>% 
    group_by(CancerType) %>% 
    summarize(start=min(id), end=max(id)) %>% 
    rowwise() %>% 
    mutate(title=mean(c(start, end)))
  GroupCol <- rbind(GroupCol, base_data2)
}

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]


##### Create long format dataframe (final_data)

data$model <- "NA"
data$y <- NA
final_data <- data

#whole
data$y <- data$diff
data$model <- "increase_whole"
final_data <- rbind(final_data, data)

data$y <- data$Rsquared
data$model <- "basal_whole"
final_data <- rbind(final_data, data)

data$y <- 1 - data$Rsquared - data$diff
data$model <- "none_whole"
final_data <- rbind(final_data, data)

#Cancer
data$y <- data$diffCancer
data$model <- "increase_cancer"
final_data <- rbind(final_data, data)

data$y <- data$RsquaredCancer
data$model <- "basal_cancer"
final_data <- rbind(final_data, data)

data$y <- 1 - data$RsquaredCancer - data$diffCancer
data$model <- "none_cancer"
final_data <- rbind(final_data, data)

#tme
data$y <- data$diffTME
data$model <- "increase_tme"
final_data <- rbind(final_data, data)

data$y <- data$RsquaredTME
data$model <- "basal_tme"
final_data <- rbind(final_data, data)

data$y <- 1 - data$RsquaredTME - data$diffTME
data$model <- "none_tme"
final_data <- rbind(final_data, data)

#order stacks
final_data$model <- factor(final_data$model, levels = c("none_tme","increase_tme", "basal_tme", 
                                                        "none_cancer","increase_cancer", "basal_cancer",
                                                        "none_whole","increase_whole", "basal_whole"))
# Make the plot
p <- ggplot(final_data, aes(x=as.factor(id), y=y, fill=model)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", width = 1) +
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0.5, xend = start, yend = 0.5), colour = "darkgrey", alpha=0.8, size=0.2 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 1.5, xend = start, yend = 1.5), colour = "darkgrey", alpha=0.8, size=0.2 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 2.5, xend = start, yend = 2.5), colour = "darkgrey", alpha=0.8, size=0.2 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(final_data$id),3), y = c(0.5, 1.5, 2.5), label = c("0.5", "0.5", "0.5") , color="darkgray", size=c(2.5,3.75, 5), angle=0, fontface="bold", hjust=c(1,1,1.2)) +
  #scale_fill_manual(values = c("yellow", "yellow", "yellow", "darkred", "darkred", "darkred", "lightblue" , "lightblue" , "lightblue")) + 
  scale_fill_manual(values = c("white", "darkred", "lightblue", "white", "darkred", "lightblue", "white" , "darkred" , "lightblue")) + 
  ylim(-1.5,3.2) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-3,4), "cm")
  ) + 
  #annotate("text", y = -0.5, x = 0, label = "R-squared gain by Neighborhood", size=4.5) + 
  coord_polar() + 
  #geom_segment(data=base_data, aes(x = start, y = mean, xend = end, yend = mean), colour = "black", alpha=1, size=0.4 , inherit.aes = FALSE )   + 
  # Add base line information
  geom_segment(data=GroupCol, aes(x = start, y = -0.05, xend = end, yend = -0.05, col=CancerType), alpha=1, size=3, inherit.aes = FALSE )+
  scale_color_manual(values = Seurat::DiscretePalette(15, palette = "polychrome") [c(1:3, 5, 6:10)]) + 
  #geom_hline(yintercept = 1) + geom_hline(yintercept = 2) + geom_hline(yintercept = 3) +
  geom_segment(data=base_data, aes(x = start, y = 1, xend = end, yend = 1), colour = "black", alpha=1, size=0.5 , inherit.aes = FALSE )  +
  geom_segment(data=base_data, aes(x = start, y = 2, xend = end, yend = 2), colour = "black", alpha=1, size=0.5 , inherit.aes = FALSE ) + 
  geom_segment(data=base_data, aes(x = start, y = 0, xend = end, yend = 0), colour = "black", alpha=1, size=0.5 , inherit.aes = FALSE ) + 
  #geom_text(data=base_data, aes(x = title, y = -0.05, label=Hallmark), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE) + 
  geom_hline(yintercept = 3)

p
ggsave("Desktop/barplot/barplot.png", bg = "white", width = 9, height = 9)  



###################
### Circular barplot ###

# library
library(tidyverse)
library(dplyr)
library(Seurat)


#merge in one table all important information
lm_mat.basal$diff <- lm_mat.1$Rsquared - lm_mat.basal$Rsquared
lm_mat.basal$RsquaredCancer <- lm_mat.basal.Cancer$Rsquared
lm_mat.basal$diffCancer <- lm_mat.1.Cancer$Rsquared - lm_mat.basal.Cancer$Rsquared
lm_mat.basal$RsquaredTME <- lm_mat.basal.TME$Rsquared
lm_mat.basal$diffTME <- lm_mat.1.TME$Rsquared - lm_mat.basal.TME$Rsquared

########### RUN from here !!!!!!!!! #########################3

lm_mat.basal <- read.table("/home/malsibai/rnaseq/expression/Neighberhood/df_models.txt")



# Define CancerType groups

lm_mat.basal$CancerType <- NA
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("Acinar", "IC", "PC1", "PC2")] <- "Prostate"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("Breast", "BreastA", "Ductal", "DuctalFFPE", "TNBCA", "M1", "M2", "M3", "M4")] <- "Breast"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("C20", "C21", "C34", "C51", "C7")] <- "Kidney"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("cHC1T", "HCC1T", "HCC2T", "HCC5D", "ICC1L")] <- "Liver"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("Colorectal", "CRC1", "CRC2", "Intestine", "Co1", "Co2", "Co3", "Co4")] <- "Colorectal"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("DU2", "DU3", "DU12", "DU8", "DU13")] <- "Bladder"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("OV4A", "OVD1", "Ovarian", "CUP295", "OVFFPE")] <- "Ovarian"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("P259_H2A2", "P264", "P270", "P288", "P306")] <- "Pancreas"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("Glioblastoma", "UKF242T", "UKF260T", "UKF269T", "UKF275T")] <- "Glioblastoma"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("P1","P3", "P4", "P5", "P6", "P7","P8")] <- "Lung"

### Major dendogram ###

#sort hallmarks
library(dendextend)
library(circlize)

df <- data.frame(matrix(nrow = 13, ncol = 46))
colnames(df) <- unique(lm_mat.basal$Sample)
rownames(df) <- paste0("H",1:13)
for (h in paste0("H",1:13)) {
  for (sample in unique(lm_mat.basal$Sample)) {
    df[h, sample] <- lm_mat.basal[lm_mat.basal$Sample==sample & lm_mat.basal$Hallmark == h, "diff"]
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

data <- lm_mat.basal[order(lm_mat.basal$CancerType),]
data$Hallmark <- factor(data$Hallmark, levels = hallmark_order)

data <- data %>%
  group_by(CancerType, Hallmark) %>%
  #mutate(Rsquared = mean(Rsquared)) %>%
  #mutate(diff = mean(diff)) %>%
  mutate(RsquaredCancer = mean(RsquaredCancer)) %>%
  mutate(diffCancer = mean(diffCancer)) %>%
  mutate(RsquaredTME = mean(RsquaredTME)) %>%
  mutate(diffTME = mean(diffTME)) %>%
  mutate(RsquaredBuffer = mean(RsquaredBuffer)) %>%
  mutate(diffBuffer = mean(diffBuffer)) %>%
  ungroup()

data[c( "id", "EstimateSlope", "EstimateStdError", "EstimateTvalue", "EstimatePvalue")] <- NULL
data <- data[duplicated(data[,-1]) == F,]

data <- data.frame(data)

data <- select(data, -c(Rsquared, diff))

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 10
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
#tmp <- data[!is.na(data),] %>% group_by(Hallmark) %>% summarise(mean = mean(diff), sd = sd(diff))
#base_data$mean <- tmp$mean[1:13]


# get locations for each cancer group to put color lines
GroupCol <- c()
for (h in paste0("H", 1:13)) {
  base_data2 <- data[data$Hallmark==h & !is.na(data$Sample),] %>%
    group_by(CancerType) %>%
    summarize(start=min(id), end=max(id)) %>%
    rowwise() %>%
    mutate(title=mean(c(start, end)))
  GroupCol <- rbind(GroupCol, base_data2)
}


#GroupCol[GroupCol$CancerType == "Prostate",]$end <- GroupCol[GroupCol$CancerType == "Prostate",]$end  +1
#GroupCol[GroupCol$CancerType == "Prostate",]$start <- GroupCol[GroupCol$CancerType == "Prostate",]$start +1

GroupCol$end <- GroupCol$end + 1

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]


##### Create long format dataframe (final_data)

data$model <- "NA"
data$y <- NA
final_data <- data

#whole
#data$y <- data$diff
#data$model <- "increase_whole"
#final_data <- rbind(final_data, data)

#data$y <- data$Rsquared
#data$model <- "basal_whole"
#final_data <- rbind(final_data, data)

#data$y <- 1 - data$Rsquared - data$diff
#data$model <- "none_whole"
#final_data <- rbind(final_data, data)

#Cancer
data$y <- data$diffCancer
data$model <- "increase_cancer"
final_data <- rbind(final_data, data)

data$y <- data$RsquaredCancer
data$model <- "basal_cancer"
final_data <- rbind(final_data, data)

data$y <- 1 - data$RsquaredCancer - data$diffCancer
data$model <- "none_cancer"
final_data <- rbind(final_data, data)

#tme
data$y <- data$diffTME
data$model <- "increase_tme"
final_data <- rbind(final_data, data)

data$y <- data$RsquaredTME
data$model <- "basal_tme"
final_data <- rbind(final_data, data)

data$y <- 1 - data$RsquaredTME - data$diffTME
data$model <- "none_tme"
final_data <- rbind(final_data, data)


#buffer
data$y <- data$diffBuffer
data$model <- "increase_buffer"
final_data <- rbind(final_data, data)

data$y <- data$RsquaredBuffer
data$model <- "basal_buffer"
final_data <- rbind(final_data, data)

data$y <- 1 - data$RsquaredBuffer - data$diffBuffer
data$model <- "none_buffer"
final_data <- rbind(final_data, data)



#order stacks
final_data$model <- factor(final_data$model, levels = c("none_tme","increase_tme", "basal_tme",
                                                        "none_buffer","increase_buffer", "basal_buffer",
                                                        "none_cancer","increase_cancer", "basal_cancer"))



# Make the plot
p <- ggplot(final_data, aes(x=as.factor(id), y=y, fill=model)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", width = 1) +
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0.5, xend = start, yend = 0.5), colour = "darkgrey", alpha=0.8, size=0.2 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 1.5, xend = start, yend = 1.5), colour = "darkgrey", alpha=0.8, size=0.2 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 2.5, xend = start, yend = 2.5), colour = "darkgrey", alpha=0.8, size=0.2 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(final_data$id),3), y = c(0.5, 1.5, 2.5), label = c("0.5", "0.5", "0.5") , color="darkgray", size=c(2.5,2.75, 3), angle=0, fontface="bold", hjust=c(1,1,1)) +
  annotate("text", x = rep(max(final_data$id),3), y = c(0.8, 1.8, 2.8), label = c("CANCER", "Buffer", "TME") , color="darkgray", size=c(2.5,2.75, 3), angle=0, fontface="bold", hjust=c(0.7,0.8,0.7)) +
  #scale_fill_manual(values = c("yellow", "yellow", "yellow", "darkred", "darkred", "darkred", "lightblue" , "lightblue" , "lightblue")) +
  scale_fill_manual(values = c("white", "darkred", "lightblue", "white", "darkred", "lightblue", "white" , "darkred" , "lightblue")) +
  ylim(-1.5,3.2) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-3,4), "cm")
  ) +
  #annotate("text", y = -0.5, x = 0, label = "R-squared gain by Neighborhood", size=4.5) +
  coord_polar() +
  #geom_segment(data=base_data, aes(x = start, y = mean, xend = end, yend = mean), colour = "black", alpha=1, size=0.4 , inherit.aes = FALSE )   +
  # Add base line informationl
  #geom_point(data=GroupCol, aes(x = start, y = -0.05, col=CancerType), alpha=1, size=0.5, inherit.aes = FALSE) +
  geom_segment(data=GroupCol, aes(x = start, y = 0, xend = start, yend = -0.2, col=CancerType), alpha=1, size=1.5, inherit.aes = FALSE )+
  scale_color_manual(values = Seurat::DiscretePalette(15, palette = "polychrome") [c(1:3, 5, 6:10, 12)]) +
  #geom_hline(yintercept = 1) + geom_hline(yintercept = 2) + geom_hline(yintercept = 3) +
  geom_segment(data=base_data, aes(x = start - 0.5, y = 1, xend = end + 0.5, yend = 1), colour = "black", alpha=1, size=0.5 , inherit.aes = FALSE )  +
  geom_segment(data=base_data, aes(x = start - 0.5, y = 2, xend = end + 0.5, yend = 2), colour = "black", alpha=1, size=0.5 , inherit.aes = FALSE ) +
  geom_segment(data=base_data, aes(x = start - 0.5, y = 0, xend = end + 0.5, yend = 0), colour = "black", alpha=1, size=0.5 , inherit.aes = FALSE ) +
  geom_text(data=base_data, aes(x = title, y = 3.2, label=Hallmark), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE) +
  geom_hline(yintercept = 3,color= "darkgrey")

p
#ggsave("Desktop/barplot/barplot3.png", bg = "white", width = 14, height = 9)  


####################################################################

data <- lm_mat.basal[order(lm_mat.basal$CancerType),]
data$Hallmark <- factor(data$Hallmark, levels = hallmark_order)

data <- data %>%
  group_by(CancerType, Hallmark) %>%
  mutate(Rsquared = mean(Rsquared)) %>%
  mutate(diff = mean(diff)) %>%
  ungroup()

data[c( "id", "EstimateSlope", "EstimateStdError", "EstimateTvalue", "EstimatePvalue")] <- NULL
data <- data[duplicated(data[,-1]) == F,]

data <- data.frame(data)

#data <- select(data, -c(Rsquared, diff))

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 10
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
#tmp <- data[!is.na(data),] %>% group_by(Hallmark) %>% summarise(mean = mean(diff), sd = sd(diff))
#base_data$mean <- tmp$mean[1:13]


# get locations for each cancer group to put color lines
GroupCol <- c()
for (h in paste0("H", 1:13)) {
  base_data2 <- data[data$Hallmark==h & !is.na(data$Sample),] %>%
    group_by(CancerType) %>%
    summarize(start=min(id), end=max(id)) %>%
    rowwise() %>%
    mutate(title=mean(c(start, end)))
  GroupCol <- rbind(GroupCol, base_data2)
}


#GroupCol[GroupCol$CancerType == "Prostate",]$end <- GroupCol[GroupCol$CancerType == "Prostate",]$end  +1
#GroupCol[GroupCol$CancerType == "Prostate",]$start <- GroupCol[GroupCol$CancerType == "Prostate",]$start +1

GroupCol$end <- GroupCol$end + 1

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]


##### Create long format dataframe (final_data)

data$model <- "NA"
data$y <- NA
final_data <- data

#whole
data$y <- data$diff
data$model <- "increase_whole"
final_data <- rbind(final_data, data)

data$y <- data$Rsquared
data$model <- "basal_whole"
final_data <- rbind(final_data, data)

data$y <- 1 - data$Rsquared - data$diff
data$model <- "none_whole"
final_data <- rbind(final_data, data)


#order stacks
final_data$model <- factor(final_data$model, levels = c("none_tme","increase_tme", "basal_tme",
                                                        "none_buffer","increase_buffer", "basal_buffer",
                                                        "none_cancer","increase_cancer", "basal_cancer"))



# Make the plot
p <- ggplot(final_data, aes(x=as.factor(id), y=y, fill=model)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", width = 1) +
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0.5, xend = start, yend = 0.5), colour = "darkgrey", alpha=0.8, size=0.2 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 1.5, xend = start, yend = 1.5), colour = "darkgrey", alpha=0.8, size=0.2 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 2.5, xend = start, yend = 2.5), colour = "darkgrey", alpha=0.8, size=0.2 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(final_data$id),3), y = c(0.5, 1.5, 2.5), label = c("0.5", "0.5", "0.5") , color="darkgray", size=c(2.5,2.75, 3), angle=0, fontface="bold", hjust=c(1,1,1)) +
  annotate("text", x = rep(max(final_data$id),3), y = c(0.8, 1.8, 2.8), label = c("Whole", "Buffer", "TME") , color="darkgray", size=c(2.5,2.75, 3), angle=0, fontface="bold", hjust=c(0.7,0.8,0.7)) +
  #scale_fill_manual(values = c("yellow", "yellow", "yellow", "darkred", "darkred", "darkred", "lightblue" , "lightblue" , "lightblue")) +
  scale_fill_manual(values = c("white", "darkred", "lightblue", "white", "darkred", "lightblue", "white" , "darkred" , "lightblue")) +
  ylim(-1.5,1.2) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-3,4), "cm")
  ) +
  #annotate("text", y = -0.5, x = 0, label = "R-squared gain by Neighborhood", size=4.5) +
  coord_polar() +
  #geom_segment(data=base_data, aes(x = start, y = mean, xend = end, yend = mean), colour = "black", alpha=1, size=0.4 , inherit.aes = FALSE )   +
  # Add base line informationl
  #geom_point(data=GroupCol, aes(x = start, y = -0.05, col=CancerType), alpha=1, size=0.5, inherit.aes = FALSE) +
  geom_segment(data=GroupCol, aes(x = start, y = 0, xend = start, yend = -0.2, col=CancerType), alpha=1, size=4, inherit.aes = FALSE )+
  scale_color_manual(values = Seurat::DiscretePalette(15, palette = "polychrome") [c(1:3, 5, 6:10, 12)]) +
  #geom_hline(yintercept = 1) + geom_hline(yintercept = 2) + geom_hline(yintercept = 3) +
  geom_segment(data=base_data, aes(x = start - 0.5, y = 0, xend = end + 0.5, yend = 0), colour = "black", alpha=1, size=2 , inherit.aes = FALSE ) +
  geom_text(data=base_data, aes(x = title, y = 3.2, label=Hallmark), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE) +
  geom_hline(yintercept = 1,color= "darkgrey")


p
ggsave("Desktop/barplot/barplotwhole.png", bg = "white", width = 14, height = 9)  
