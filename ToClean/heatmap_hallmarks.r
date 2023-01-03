
ST_enhanced.St <- readRDS("Desktop/IJC/datasets/Public/IC/RDS/")
data_name <- "IC"
ST_enhanced.St <- readRDS(paste0("Desktop/IJC/datasets/Public/", data_name, "/RDS/", data_name, "_enhanced.rds"))

colnames(ST_enhanced.St@meta.data)
H_modules <- ST_enhanced.St@meta.data[,18:31]
#H_modules <- ST_enhanced.St@meta.data[,c(44, 13:25)]
H_modules[2:14] <- scale(H_modules[2:14])
mean(H_modules$H12)
sd(H_modules$H12)

H_modules$NMF_clusters <-  paste0(data_name, "_", H_modules$NMF_clusters)
H_modules$Cancer <- data_name

rownames(H_modules) <- paste0(data_name, "_", rownames(H_modules)) 
H_modules$subspot <- rownames(H_modules)
head(H_modules)


write.table(H_modules, paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hallmarks/", data_name, ".txt"), sep = "\t")

a <- read.table("Desktop/IJC/datasets/IGTP/figuresPaper/hallmarks/Ovarian.txt", sep = "\t")





files <- list.files(path="Desktop/IJC/datasets/IGTP/figuresPaper/hallmarks", pattern="*.txt", full.names=TRUE, recursive=FALSE)

all <- lapply(files, function(x){
  mat <- read.table(x, sep = "\t")
})

library(rlist)
H_modules <- list.rbind(all)

matrix <- read.table("Desktop/H_modules_Liver_Acinar.xls", sep = "\t", header = T)
rownames(matrix) <- matrix$subspot

H_modules <- rbind(matrix, H_modules)
unique(H_modules$NMF_clusters)
H_modules$TME_specific <- "a"

H_modules$TME_specific[H_modules$NMF_clusters == "Breast_1"] <- "ECM"
H_modules$TME_specific[H_modules$NMF_clusters == "Breast_2"] <- "ECM"
H_modules$TME_specific[H_modules$NMF_clusters == "Breast_3"] <- "Mix"
H_modules$TME_specific[H_modules$NMF_clusters == "Breast_4"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "Breast_5"] <- "Tumor"

H_modules$TME_specific[H_modules$NMF_clusters == "OV4A_1"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "OV4A_2"] <- "Immune"
H_modules$TME_specific[H_modules$NMF_clusters == "OV4A_3"] <- "ECM"
H_modules$TME_specific[H_modules$NMF_clusters == "OV4A_4"] <- "Endothelial"
H_modules$TME_specific[H_modules$NMF_clusters == "OV4A_5"] <- "ECM"

H_modules$TME_specific[H_modules$NMF_clusters == "NASH_1"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "NASH_2"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "NASH_3"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "NASH_4"] <- "ECM"
H_modules$TME_specific[H_modules$NMF_clusters == "NASH_5"] <- "Tumor"

H_modules$TME_specific[H_modules$NMF_clusters == "HBV_1"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "HBV_2"] <- "Immune"
H_modules$TME_specific[H_modules$NMF_clusters == "HBV_3"] <- "Unknown"
H_modules$TME_specific[H_modules$NMF_clusters == "HBV_4"] <- "Unknown"
H_modules$TME_specific[H_modules$NMF_clusters == "HBV_5"] <- "Tumor"

H_modules$TME_specific[H_modules$NMF_clusters == "HCV1_1"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "HCV1_2"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "HCV1_3"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "HCV1_4"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "HCV1_5"] <- "Immune"

H_modules$TME_specific[H_modules$NMF_clusters == "Ductal_1"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "Ductal_2"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "Ductal_3"] <- "Immune"
H_modules$TME_specific[H_modules$NMF_clusters == "Ductal_4"] <- "Unknown"
H_modules$TME_specific[H_modules$NMF_clusters == "Ductal_5"] <- "Immune"

H_modules$TME_specific[H_modules$NMF_clusters == "HCV2_1"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "HCV2_2"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "HCV2_3"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "HCV2_4"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "HCV2_5"] <- "Tumor"

H_modules$TME_specific[H_modules$NMF_clusters == "Ovarian_1"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "Ovarian_2"] <- "ECM"
H_modules$TME_specific[H_modules$NMF_clusters == "Ovarian_3"] <- "ECM"
H_modules$TME_specific[H_modules$NMF_clusters == "Ovarian_4"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "Ovarian_5"] <- "Immune"

H_modules$TME_specific[H_modules$NMF_clusters == "Glioblastoma_1"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "Glioblastoma_2"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "Glioblastoma_3"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "Glioblastoma_4"] <- "ECM"
H_modules$TME_specific[H_modules$NMF_clusters == "Glioblastoma_5"] <- "Tumor"

H_modules$TME_specific[H_modules$NMF_clusters == "Acinar_1"] <- "ECM"
H_modules$TME_specific[H_modules$NMF_clusters == "Acinar_2"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "Acinar_3"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "Acinar_4"] <- "Immune"
H_modules$TME_specific[H_modules$NMF_clusters == "Acinar_5"] <- "Tumor"

H_modules$TME_specific[H_modules$NMF_clusters == "IC_1"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "IC_2"] <- "ECM"
H_modules$TME_specific[H_modules$NMF_clusters == "IC_3"] <- "Muscle"
H_modules$TME_specific[H_modules$NMF_clusters == "IC_4"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "IC_5"] <- "Tumor"

H_modules$TME_specific[H_modules$NMF_clusters == "Colorectal_1"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "Colorectal_2"] <- "ECM"
H_modules$TME_specific[H_modules$NMF_clusters == "Colorectal_3"] <- "Mix"
H_modules$TME_specific[H_modules$NMF_clusters == "Colorectal_4"] <- "Tumor"
H_modules$TME_specific[H_modules$NMF_clusters == "Colorectal_5"] <- "Immune"


H_modules$TME_general <- H_modules$TME_specific
H_modules$TME_general[!H_modules$TME_specific %in% c("Tumor", "Immune")] <- "Stroma"


#sample 
H_modules_sampled <- c()

for (cancer_type in unique(H_modules$Cancer)) {
  cancer_sub <- H_modules[H_modules$Cancer == cancer_type,]
  H_modules_sampled <- rbind(H_modules_sampled, cancer_sub[sample(rownames(cancer_sub), size = round(nrow(cancer_sub)*0.2), replace = F),])
}


#### Heatmap
library(pheatmap)
library(RColorBrewer)
library(dplyr)


H_modules_avg <- arrange(H_modules_sampled, Cancer)

ann_heatmap <- H_modules_avg[c("NMF_clusters", "TME_specific", "TME_general","Cancer")]
ann_heatmap <- data.frame(ann_heatmap, row.names = 1)
#ann_heatmap <- select(ann_heatmap, "Cancer")

H_modules_avg.mt <- H_modules_avg[,c(1,3:15)]
H_modules_avg.mt <- data.matrix(data.frame(H_modules_avg.mt, row.names = 1))

palette <- RColorBrewer::brewer.pal(8, name = "Paired")
names(palette) <- unique(H_modules_avg$Cancer)
ann_colors <- list(TME_specific = c(Tumor="darkred", Immune="darkgreen", ECM="darkorange1", Unknown="gray49", Mix="darkorchid3", Muscle="red", Endothelial="yellow"),
                   TME_general = c(Tumor="darkred", Immune="darkgreen", Stroma="darkblue"),
                Cancer = palette)


pheatmap(t(H_modules_avg.mt),
         scale = "none",
         annotation_col = ann_heatmap,
         annotation_colors = ann_colors,
         legend = T,
         cluster_cols = T,
         cluster_rows = T,
         show_rownames = T,
         show_colnames = F, clustering_distance_cols = "manhattan", clustering_distance_rows = "manhattan",
         fontsize_col = 2)

H_modules_avg <- H_modules %>% group_by(NMF_clusters, Cancer) %>% summarise(across(H1:H13, ~ mean(.x)))

                                                                    



H_modules_avg$TME_specific <- "a"

H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Breast_1"] <- "ECM"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Breast_2"] <- "ECM"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Breast_3"] <- "Mix"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Breast_4"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Breast_5"] <- "Tumor"

H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "OV4A_1"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "OV4A_2"] <- "Immune"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "OV4A_3"] <- "ECM"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "OV4A_4"] <- "Endothelial"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "OV4A_5"] <- "ECM"

H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "NASH_1"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "NASH_2"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "NASH_3"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "NASH_4"] <- "ECM"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "NASH_5"] <- "Tumor"

H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "HBV_1"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "HBV_2"] <- "Immune"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "HBV_3"] <- "Unknown"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "HBV_4"] <- "Unknown"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "HBV_5"] <- "Tumor"

H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "HCV1_1"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "HCV1_2"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "HCV1_3"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "HCV1_4"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "HCV1_5"] <- "Immune"

H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Ductal_1"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Ductal_2"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Ductal_3"] <- "Immune"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Ductal_4"] <- "Unknown"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Ductal_5"] <- "Immune"

H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "HCV2_1"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "HCV2_2"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "HCV2_3"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "HCV2_4"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "HCV2_5"] <- "Tumor"

H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Ovarian_1"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Ovarian_2"] <- "ECM"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Ovarian_3"] <- "ECM"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Ovarian_4"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Ovarian_5"] <- "Immune"

H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Glioblastoma_1"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Glioblastoma_2"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Glioblastoma_3"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Glioblastoma_4"] <- "ECM"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Glioblastoma_5"] <- "Tumor"

H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Acinar_1"] <- "ECM"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Acinar_2"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Acinar_3"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Acinar_4"] <- "Immune"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Acinar_5"] <- "Tumor"

H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "IC_1"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "IC_2"] <- "ECM"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "IC_3"] <- "Muscle"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "IC_4"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "IC_5"] <- "Tumor"

H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Colorectal_1"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Colorectal_2"] <- "ECM"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Colorectal_3"] <- "Mix"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Colorectal_4"] <- "Tumor"
H_modules_avg$TME_specific[H_modules_avg$NMF_clusters == "Colorectal_5"] <- "Immune"


H_modules_avg$TME_general <- H_modules_avg$TME_specific
H_modules_avg$TME_general[!H_modules_avg$TME_specific %in% c("Tumor", "Immune")] <- "Stroma"

