##################################################
## Project: Cancer Hallmarks
## Script purpose: Plot a heatmap representing and the hallmark activity within each ESTIMATE cluster
## Author: Sergi Cervilla* & Mustafa Sibai*
##################################################

library(pheatmap)
library(paletteer)
library(estimate)
library(dplyr)
library(Seurat)
library(ggpubr)
library(ClassDiscovery)



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



# Significance test of the obtained clustering of the cancer type based on scaled hallmarks activities
distance_mat <- ClassDiscovery::distanceMatrix(t(H_EST.mt), metric = "pearson")
Hierar_cl <- hclust(distance_mat, method = "average")

fit <- cutree(Hierar_cl, k = 2)

t <- table(sapply(spot_info$sample, function(x) {
  rep(annotation_tissue[[x]], )
}), fit)[-4,]
chisq.test(t)


