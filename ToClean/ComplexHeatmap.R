mt_importances <- read.table("Downloads/RF_Imp.mt.hmap.cancer.txt")
head(mt_importances)

#create same clusters as in the heatmap
distance_mat <- ClassDiscovery::distanceMatrix(t(mt_importances), metric = "maximum")
Hierar_cl <- hclust(distance_mat, method = "average")

plot(Hierar_cl)
fit <- cutree(Hierar_cl, k = 8)
mt_importances$fit <- fit
rect.hclust(Hierar_cl, k=8, border="green")


#Add metadata to the matrix: Tissue type and sample
mt_importances$tissue <- sapply(rownames(mt_importances), function(x) {
  if (str_split(x, pattern = "_")[[1]][1] == "P259") return("Pancreas")
  return(tissue[[str_split(x, pattern = "_")[[1]][1]]])
})

mt_importances$sample <- sapply(rownames(mt_importances), function(x) {
  if (str_split(x, pattern = "_")[[1]][1] == "P259") return("P259_H2A2")
  return(str_split(x, pattern = "_")[[1]][1])
})

annot <- mt_importances[,c("fit", "tissue")]
annot$fit <- factor(annot$fit)
pheatmap(t(mt_importances[,1:7]), clustering_distance_cols = "maximum", color = rev(paletteer_c("grDevices::Reds 2", 30) ), cutree_cols = 8,
         clustering_method = "average", show_colnames = F, cluster_rows = F, annotation_col = annot)
pheatmap(t/colSums(t))


t <- table(mt_importances$tissue, mt_importances$fit)
t <- table(mt_importances$sample, mt_importances$fit)
pheatmap(t/colSums(t))
fisher.test(t )
chisq.test(t[,1:2])

###
#count the percentage of tissue type in each cluster
counts <- lapply(1:8,function(fit) {
 sub <-  mt_importances[mt_importances$fit==fit, c("tissue", "sample")]
 #size <- nrow(sub)
 table(factor(sapply(unique(sub$sample), function(x) {tissue[[x]]}),levels=unique(mt_importances$tissue)))/length(unique(sub$sample))
})
#convert list to data frame
counts <- do.call(cbind, counts)
colnames(counts) <- 1:8

#add count values for each of model (columns)
ann_counts_table <- mt_importances %>% select(fit)
for (i in colnames(counts)) {
  for(cancer in rownames(counts)) {
    ann_counts_table[ann_counts_table$fit == i, cancer] <- counts[cancer, i]
  }
}

#Creation of the heatmap
library(ComplexHeatmap)
#color mapping for 
imp_col = circlize::colorRamp2(c(-1, 0, 1), c("blue","white", "red")) 
#annotation on the top area
ann_top <-  HeatmapAnnotation(
    H1  = anno_simple(rnorm(290), col=imp_col), #first argument should be the correlations of predictor/SHAPLEE
    H13  = anno_simple(rnorm(290), col=imp_col),
    H3  = anno_simple(rnorm(290), col=imp_col),
    H5  = anno_simple(rnorm(290), col=imp_col),
    H6  = anno_simple(rnorm(290), col=imp_col),
    H7  = anno_simple(rnorm(290), col=imp_col),
    H9  = anno_simple(rnorm(290), col=imp_col),
  annotation_name_gp= gpar(fontsize = 10), #size of labels
  annotation_name_side = "right", #position of labels
  height = unit(1.2, "cm") #size of the annotation layer
  )

#color map for tissue percentage
perc_col = circlize::colorRamp2(c(0, 1), c("white", "black")) 
ann <-  HeatmapAnnotation(
  Bladder  = anno_simple(ann_counts_table$Bladder, col=perc_col),
  Breast  = anno_simple(ann_counts_table$Breast, col=perc_col),
  Colorectal  = anno_simple(ann_counts_table$Colorectal, col=perc_col),
  Glioblastoma  = anno_simple(ann_counts_table$Glioblastoma, col=perc_col),
  kidney  = anno_simple(ann_counts_table$Kidney, col=perc_col),
  Liver  = anno_simple(ann_counts_table$Liver, col=perc_col),
  Lung  = anno_simple(ann_counts_table$Lung, col=perc_col),
  Ovarian  = anno_simple(ann_counts_table$Ovarian, col=perc_col),
  Pancreas  = anno_simple(ann_counts_table$Pancreas, col=perc_col),
  Prostate  = anno_simple(ann_counts_table$Prostate, col=perc_col),
  annotation_name_side = "right", #position of labels
  height = unit(3, "cm"), #size of the annotation layer
  annotation_name_gp= gpar(fontsize = 7), #size of labels
)



Heatmap(t(mt_importances[,1:7]), 
        col=rev(paletteer_c("grDevices::Reds 2", 30) ), #colors
        cluster_rows = F, clustering_distance_columns = "maximum", show_column_names = F, #clustering (similar to pheatmap)
        clustering_method_columns = "average", 
        top_annotation = ann_top, bottom_annotation = ann, #annotations
        column_split = 8) #dendrogram split
#add a bottom line to separate annotation (I haven't explored this yet)
decorate_annotation("Bladder", {
  grid.lines(unit(c(-40, 400), "mm"), unit(c(1, 1), "npc"))
})
