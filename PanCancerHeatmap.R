##################################################
## Project: Cancer Hallmarks
## Script purpose: Plot a heatmap representing and the hallmark activity within each ESTIMATE cluster
## Date: 22/12/2022
## Author: Sergi Cervilla & Mustafa Sibai
##################################################

library(pheatmap)
library(paletteer)
library(estimate)
library(dplyr)



### Pan-Cancer analysis using Estimate clusters

#Take name of the samples (for instance in ESTIMATE clusters files)
files <- grep(list.files("..", full.names = F,), pattern='*.sh', invert=TRUE, value=TRUE)
samples <- c()
for (file in files) {
  x <- strsplit(file, split = ".txt")[[1]][1]
  samples <- c(samples, x)
}

#load average activity for each cluster and sample
H_EST <- lapply(samples, function(sample) {
  tmp <- read.table("", sep = "\t", header = T)
  tmp$tumor <- sample
  return(tmp)
})
H_EST <- do.call(rbind, H_EST)
#add cluster label
H_EST$cluster <- factor(rep(1:5, 58))

source("../utils/SamplesMetadata.R")

#unique row ID
H_EST$ID <- rownames(H_EST)
#add used embedding method metadata
H_EST$embedding_method <- ""
for (n in unique(H_EST$tumor)) {
  H_EST$embedding_method[H_EST$tumor == n] <- annotation_method[n]
}
H_EST$embedding_method <- as.character(H_EST$embedding_method)

#add anatomic tumor site metadata
H_EST$Anatomic_site <- ""
for (n in unique(H_EST$tumor)) {
  H_EST$Anatomic_site[H_EST$tumor==n] <- annotation_site[n]
}
H_EST$Anatomic_site <- as.character(H_EST$Anatomic_site)


ann_heatmap <- H_EST[c( "ID", "embedding_method", "Anatomic_site", "cluster")]
ann_heatmap <- data.frame(ann_heatmap, row.names = 1)


#Extract ID and average Hallmark activities
H_EST.mt <- H_EST[,c(16,1:13)]
#Put ID as rownames
H_EST.mt <- data.matrix(data.frame(H_EST.mt, row.names = 1))
#Use full hallmark names in the matrix
colnames(H_EST.mt) <- hallmark_names

#Color for each anatomic tumor site
palette <- Seurat::DiscretePalette(15, palette = "polychrome") [c(1:3, 5, 6:10, 12)]
names(palette) <- sort(unique(ann_heatmap$Anatomic_site))
#Color for FFPE samples
FFPE <- Seurat::DiscretePalette(15, palette = "polychrome")[11]
#Color for OCT samples
OCT <- Seurat::DiscretePalette(15, palette = "polychrome")[14]
#Color for ESTIMATE clusters
cluster <- c("5" = "lightgoldenrod1", "4" = "lightgoldenrod3", "3" = "lightpink2", "2"= "orchid3", "1"= "orchid4")
#colors for the annotation layer
ann_colors <- list(embedding_method = c(FFPE=FFPE, OCT=OCT), Anatomic_site=palette, cluster=cluster)

#plot Pan-Cancer heatmap
pheatmap(t(H_EST.mt),
         scale = "none",
         annotation_col = ann_heatmap,
         annotation_colors = ann_colors,
         legend = T,
         cluster_cols = T,
         cluster_rows = T,
         show_rownames = T,
         show_colnames = F, clustering_distance_cols = "correlation", clustering_distance_rows = "correlation", clustering_method = "average",
         fontsize_row = 15, border_color = "black", color =   rev(paletteer_c("ggthemes::Orange-Blue Diverging", 50) ),
         cutree_cols = 2, cutree_rows = 2)



################################################################################

#6.5 x 18 inches
library(ggpubr)

ggboxplot(H_EST, x= "cluster", y = "H12") +stat_compare_means()
ggboxplot(H_EST, x= "cluster", y = "H11", facet.by = "Anatomic_site", scales = "free")
ggboxplot(H_EST, x= "cluster", y = "H13", facet.by = "embedding_method", scales = "free")+stat_compare_means()





################################################################################
annotation_compartments <- sapply(rownames(ann_heatmap), function(cluster){
  if (ann_heatmap[cluster, "fit"] == 3 & ann_heatmap[cluster, "Kmeans"] == 2)  return("Cancer")
  else if (ann_heatmap[cluster, "fit"] == 2 & ann_heatmap[cluster, "Kmeans"] == 1) return("TME")
  else (return("Buffer"))
})

saveRDS(annotation_compartments, "")
annotation_compartments <- readRDS("")

#### boxplots

data <- t(H_EST.mt)

distance_mat <- ClassDiscovery::distanceMatrix(data, metric = "pearson")
Hierar_cl <- hclust(distance_mat, method = "average")

plot(Hierar_cl)

fit <- cutree(Hierar_cl, k = 3)
fit
rect.hclust(Hierar_cl, k=3, border="green")



H_EST$fit <- fit

H_EST$scores <- df_BS[rownames(H_EST),3]

ggboxplot(H_EST, x = "fit", y = "scores")


write.table(H_EST, "", quote = F, row.names = F, sep = "\t")

######################## max and min activity per sample for each hallmark ###############
H_stats <- data.frame(tumor = rep(unique(H_EST$tumor), each = 2),
                      stat = bind_rows(rep(list(tibble(stat = c("Max", "Min"))), 41)))

for (hallmark in paste0("H", 1:13)) {
  H_stats[,hallmark] <- ""
}

for(sample in unique(H_EST$tumor)) {
  table <- filter(H_EST, tumor == sample)
  for (hallmark in paste0("H", 1:13)){
    for(c in 1:nrow(table)){
      if (table[c,hallmark] == max(table[,hallmark])) {
        H_stats[H_stats$tumor == sample & H_stats$stat == "Max",][hallmark] <- table[c, "cluster"]
      } else if (table[c,hallmark] == min(table[,hallmark])){
        H_stats[H_stats$tumor == sample & H_stats$stat == "Min",][hallmark] <- table[c, "cluster"]
      }
    }
  }
}

barplot(table(as.numeric(H_stats[H_stats$stat=="Max",]$H1))) +title("H1")
barplot(table(as.numeric(H_stats[H_stats$stat=="Max",]$H2))) +title("H2")
barplot(table(as.numeric(H_stats[H_stats$stat=="Max",]$H3))) +title("H3")
barplot(table(as.numeric(H_stats[H_stats$stat=="Max",]$H4))) +title("H4")
barplot(table(as.numeric(H_stats[H_stats$stat=="Max",]$H5))) +title("H5")
barplot(table(as.numeric(H_stats[H_stats$stat=="Max",]$H6))) +title("H6")
barplot(table(as.numeric(H_stats[H_stats$stat=="Max",]$H7))) +title("H7")
barplot(table(as.numeric(H_stats[H_stats$stat=="Max",]$H8))) +title("H8")
barplot(table(as.numeric(H_stats[H_stats$stat=="Max",]$H9))) +title("H9")
barplot(table(as.numeric(H_stats[H_stats$stat=="Max",]$H10))) +title("H10")
barplot(table(as.numeric(H_stats[H_stats$stat=="Max",]$H11))) +title("H11")
barplot(table(as.numeric(H_stats[H_stats$stat=="Max",]$H12))) +title("H12")
barplot(table(as.numeric(H_stats[H_stats$stat=="Max",]$H13))) +title("H13")

H1.samples.max <- H_stats[H_stats$stat == "Max" & H_stats$H1 == "5",]$tumor # in 5
H2.samples.max <- H_stats[H_stats$stat == "Max" & H_stats$H2 == "1",]$tumor # in 1
H3.samples.max <- H_stats[H_stats$stat == "Max" & H_stats$H3 == "5",]$tumor # in 5
H4.samples.max <- H_stats[H_stats$stat == "Max" & H_stats$H4 == "1",]$tumor # in 1
H5.samples.max <- H_stats[H_stats$stat == "Max" & H_stats$H5 == "5",]$tumor # in 5
H6.samples.max <- H_stats[H_stats$stat == "Max" & H_stats$H6 == "5",]$tumor # in 5
H7.samples.max <- H_stats[H_stats$stat == "Max" & H_stats$H7 == "5",]$tumor # in 5
H8.samples.max <- H_stats[H_stats$stat == "Max" & H_stats$H8 == "1",]$tumor # in 1
H9.samples.max <- H_stats[H_stats$stat == "Max" & H_stats$H9 == "5",]$tumor # in 5
H10.samples.max.1 <- H_stats[H_stats$stat == "Max" & H_stats$H10 == "1",]$tumor # in 1
H10.samples.max.4 <- H_stats[H_stats$stat == "Max" & H_stats$H10 == "4",]$tumor # in 4
H11.samples.max.1 <- H_stats[H_stats$stat == "Max" & H_stats$H11 == "1",]$tumor # in 1
H11.samples.max.5 <- H_stats[H_stats$stat == "Max" & H_stats$H11 == "5",]$tumor # in 5
H12.samples.max <- H_stats[H_stats$stat == "Max" & H_stats$H12 == "1",]$tumor # in 1
H13.samples.max <- H_stats[H_stats$stat == "Max" & H_stats$H13 == "5",]$tumor # in 5