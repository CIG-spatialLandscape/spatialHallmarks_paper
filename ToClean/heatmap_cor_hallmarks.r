

#extract all the Seurat objects from the folder
files <- list.files(path="Desktop/IJC/datasets/Public/EnhancedObjects", pattern="*.rds", full.names=TRUE, recursive=FALSE)

H_cor <- c()
#iterate for each file
for (file in files) {
  print(file) #print to show in which step is executing
  STobject <- readRDS(file) #load the file
  embeddings <- STobject@reductions$NMF@cell.embeddings #extract factors embeddings
  hallmarks <- STobject@meta.data[, paste0("H", 11)] #extract hallmarks scores
  mat <- matrix(0, 5, 13) #create a matrix of 0s; dim (factors*hallmarks)
  #name columns and rownames
  colnames(mat) <- colnames(hallmarks)
  rownames(mat) <- paste0(STobject@images[[1]]@key, "_", 1:5)
  #fill each position of the matrix
  for (h in 1:13) {
    for (factor in 1:5) {
      mat[factor, h] <- cor(embeddings[,factor], hallmarks[,h]) #compute the pearson correlation between factor activity and hallmarks
    }
  }
  
  mat <- as.data.frame(mat) #create a data frame from the matrix
  mat$Cancer <- STobject@images[[1]]@key #add cancer type of the metadata
  H_cor <- rbind(H_cor, mat) #unify with the rest of correlations
  rm(STobject) #clean memory
  gc()
}





H_cor <- read.table("Desktop/IJC/datasets/IGTP/figuresPaper/CorrelationHallmarks.txt", header = T)
H_cor$NMF_clusters <- rownames(H_cor)

#Label each NMF cluster by its corresponding group
H_cor$Clusters[H_cor$NMF_clusters == "Breast_1"] <- "ECM"
H_cor$Clusters[H_cor$NMF_clusters == "Breast_2"] <- "ECM"
H_cor$Clusters[H_cor$NMF_clusters == "Breast_3"] <- "ECM-Immune"
H_cor$Clusters[H_cor$NMF_clusters == "Breast_4"] <- "Cancer"
H_cor$Clusters[H_cor$NMF_clusters == "Breast_5"] <- "Cancer"

H_cor$Clusters[H_cor$NMF_clusters == "OV_4A_1"] <- "Cancer"
H_cor$Clusters[H_cor$NMF_clusters == "OV_4A_2"] <- "Immune"
H_cor$Clusters[H_cor$NMF_clusters == "OV_4A_3"] <- "ECM"
H_cor$Clusters[H_cor$NMF_clusters == "OV_4A_4"] <- "ECM"
H_cor$Clusters[H_cor$NMF_clusters == "OV_4A_5"] <- "ECM"

H_cor$Clusters[H_cor$NMF_clusters == "NASH_1"] <- "Cancer"
H_cor$Clusters[H_cor$NMF_clusters == "NASH_2"] <- "Cancer"
H_cor$Clusters[H_cor$NMF_clusters == "NASH_3"] <- "ECM-Immune"
H_cor$Clusters[H_cor$NMF_clusters == "NASH_4"] <- "Infiltrated"
H_cor$Clusters[H_cor$NMF_clusters == "NASH_5"] <- "Infiltrated"

H_cor$Clusters[H_cor$NMF_clusters == "HBV_1"] <- "Cancer"
H_cor$Clusters[H_cor$NMF_clusters == "HBV_2"] <- "Immune"
H_cor$Clusters[H_cor$NMF_clusters == "HBV_3"] <- "ECM-Immune"
H_cor$Clusters[H_cor$NMF_clusters == "HBV_4"] <- "Infiltrated"
H_cor$Clusters[H_cor$NMF_clusters == "HBV_5"] <- "Infiltrated"

H_cor$Clusters[H_cor$NMF_clusters == "HCV1_1"] <- "Cancer"
H_cor$Clusters[H_cor$NMF_clusters == "HCV1_2"] <- "Cancer"
H_cor$Clusters[H_cor$NMF_clusters == "HCV1_3"] <- "Cancer"
H_cor$Clusters[H_cor$NMF_clusters == "HCV1_4"] <- "Infiltrated"
H_cor$Clusters[H_cor$NMF_clusters == "HCV1_5"] <- "Infiltrated"

H_cor$Clusters[H_cor$NMF_clusters == "Ductal_1"] <- "Cancer"
H_cor$Clusters[H_cor$NMF_clusters == "Ductal_2"] <- "Cancer"
H_cor$Clusters[H_cor$NMF_clusters == "Ductal_3"] <- "Immune"
H_cor$Clusters[H_cor$NMF_clusters == "Ductal_4"] <- "Cancer"
H_cor$Clusters[H_cor$NMF_clusters == "Ductal_5"] <- "Immune"

H_cor$Clusters[H_cor$NMF_clusters == "HCV2_1"] <- "Cancer"
H_cor$Clusters[H_cor$NMF_clusters == "HCV2_2"] <- "ECM"
H_cor$Clusters[H_cor$NMF_clusters == "HCV2_3"] <- "Infiltrated"
H_cor$Clusters[H_cor$NMF_clusters == "HCV2_4"] <- "Infiltrated"
H_cor$Clusters[H_cor$NMF_clusters == "HCV2_5"] <- "Infiltrated"

H_cor$Clusters[H_cor$NMF_clusters == "Ovarian_1"] <- "Cancer"
H_cor$Clusters[H_cor$NMF_clusters == "Ovarian_2"] <- "ECM-Immune"
H_cor$Clusters[H_cor$NMF_clusters == "Ovarian_3"] <- "ECM"
H_cor$Clusters[H_cor$NMF_clusters == "Ovarian_4"] <- "Cancer"
H_cor$Clusters[H_cor$NMF_clusters == "Ovarian_5"] <- "ECM-Immune"

H_cor$Clusters[H_cor$NMF_clusters == "Glioblastoma_1"] <- "Cancer"
H_cor$Clusters[H_cor$NMF_clusters == "Glioblastoma_2"] <- "ECM"
H_cor$Clusters[H_cor$NMF_clusters == "Glioblastoma_3"] <- "Cancer"
H_cor$Clusters[H_cor$NMF_clusters == "Glioblastoma_4"] <- "ECM"
H_cor$Clusters[H_cor$NMF_clusters == "Glioblastoma_5"] <- "Cancer"

H_cor$Clusters[H_cor$NMF_clusters == "Acinar_1"] <- "Muscle"
H_cor$Clusters[H_cor$NMF_clusters == "Acinar_2"] <- "Cancer"
H_cor$Clusters[H_cor$NMF_clusters == "Acinar_3"] <- "Cancer"
H_cor$Clusters[H_cor$NMF_clusters == "Acinar_4"] <- "Immune"
H_cor$Clusters[H_cor$NMF_clusters == "Acinar_5"] <- "Infiltrated"

H_cor$Clusters[H_cor$NMF_clusters == "IC_1"] <- "Cancer"
H_cor$Clusters[H_cor$NMF_clusters == "IC_2"] <- "ECM-Immune"
H_cor$Clusters[H_cor$NMF_clusters == "IC_3"] <- "Muscle"
H_cor$Clusters[H_cor$NMF_clusters == "IC_4"] <- "Muscle"
H_cor$Clusters[H_cor$NMF_clusters == "IC_5"] <- "Cancer"

H_cor$Clusters[H_cor$NMF_clusters == "Colorectal_1"] <- "Cancer"
H_cor$Clusters[H_cor$NMF_clusters == "Colorectal_2"] <- "ECM"
H_cor$Clusters[H_cor$NMF_clusters == "Colorectal_3"] <- "ECM-Immune"
H_cor$Clusters[H_cor$NMF_clusters == "Colorectal_4"] <- "Cancer"
H_cor$Clusters[H_cor$NMF_clusters == "Colorectal_5"] <- "Immune"

H_cor$TME_general <- H_cor$Clusters
H_cor$TME_general[!H_cor$Clusters %in% c("Cancer", "Immune", "Infiltrated")] <- "Stroma"




library(pheatmap)
library(RColorBrewer)
library(dplyr)


H_cor_filtered <- H_cor
H_cor_filtered <- H_cor[apply(H_cor[,1:13],1, max) >= 0.5,]
H_cor_filtered <- filter(H_cor_filtered, !Cancer %in% c("HCV2", "HBV", "NASH", "HCV1"))

H_cor_filtered$Cancer[H_cor_filtered$Cancer=="Breast"] <- "Invasive Lobular Carcinoma (Breast)"
H_cor_filtered$Cancer[H_cor_filtered$Cancer=="Ductal"] <- "Invasive Ductal Carcinoma (Breast)"
H_cor_filtered$Cancer[H_cor_filtered$Cancer=="IC"] <- "Adenocarcinoma (Prostate)"
H_cor_filtered$Cancer[H_cor_filtered$Cancer=="OV_4A"] <- "HGSOC (Ovary)"
H_cor_filtered$Cancer[H_cor_filtered$Cancer=="Ovarian"] <- "Endometrial Adenocarcinoma (Ovary)"
H_cor_filtered$Cancer[H_cor_filtered$Cancer=="Acinar"] <- "Acinar Cell Carcinoma (Prostate)"
H_cor_filtered$Cancer[H_cor_filtered$Cancer=="Colorectal"] <- "Invasive Adenocarcinoma (Colorectal)"

colnames(H_cor_filtered)[1:13] <- hallmark_names
colnames(H_cor_filtered)[14] <- "Tumor type"

ann_heatmap <- H_cor_filtered[c("NMF_clusters", "Clusters","Tumor type")]
ann_heatmap <- data.frame(ann_heatmap, row.names = 1)
#ann_heatmap <- select(ann_heatmap, "Cancer")

H_cor_filtered$`Tumor type` <- H_cor_filtered$Cancer 

H_cor_filtered.mt <- H_cor_filtered[,c(14,1:13)]
H_cor_filtered.mt <- data.matrix(data.frame(H_cor_filtered.mt, row.names = 1))
colnames(H_cor_filtered.mt) <- hallmark_names

palette <- RColorBrewer::brewer.pal(8, name = "Set2")
palette[2] <- RColorBrewer::brewer.pal(8, name = "Set3")[8]

names(palette) <- unique(H_cor_filtered$`Tumor type`)
ann_colors <- list(Clusters = c(Cancer="#E31A1C", Immune="#FF7F00", ECM="#33A02C",  `ECM-Immune`="#FFFF99", Muscle="#A6CEE3"),
                   `Tumor type` = palette)


pheatmap(t(H_cor_filtered.mt),
         scale = "none",
         annotation_col = ann_heatmap,
         annotation_colors = ann_colors,
         legend = T,
         cluster_cols = T,
         cluster_rows = T,
         show_rownames = T,
         show_colnames = F, clustering_distance_cols = "correlation", clustering_distance_rows = "correlation", clustering_method = "complete",
         fontsize_row = 15, border_color = "black", color = viridis::inferno(100))





