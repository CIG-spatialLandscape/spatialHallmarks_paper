##################################################
## Project: Cancer Hallmarks
## Script purpose: Generate CNV clusters and run CNV experiment
## Author: Sergi Cervilla* & Mustafa Sibai*
##################################################

library(tidyverse)
library(stringr)
library(Seurat)
library(ComplexHeatmap)
library(circlize)

#extract all samples names from files
files <- list.files("", full.names = F)
files <- stringr::str_remove(files, pattern = ".txt")
#Extract count data from the filtered sample (Spot Level)
for (sample in files) {
  #load Seurat object at spot level resolution
  STobject <- readRDS("")
  #obtain expression matrix
  mat <- as.data.frame(STobject@assays[[1]]@counts)
  #replace spot name by its coordinates (ROWxCOLUMN format)
  colnames(mat) <- paste0(STobject@images[[1]]@coordinates[,"row"], "x",STobject@images[[1]]@coordinates[,"col"])
  #create a folder for the CNV output
  dir.create(paste0("", sample))
  #save gene expression matrix
  write.table(mat, "", quote = F, sep = "\t")
  gc()
}

############################
#Run STARCH (STARCH_command.sh)
############################

#Add CNV clusters to enhanced Visium sub-spots
for (sample in files) {
  print(sample)
  #load labels of CNV clusters
  cnv_clusters <- read.csv("", row.names = 1)
  #single cell experiment enhanced object without imputed genes
  sce <- readRDS("")
  meta <- as.data.frame(colData(sce))  
  for (spot.id in rownames(cnv_clusters)) {
    coord <- as.integer(strsplit(spot.id, split = "x")[[1]])
    spots <- rownames(meta[meta$spot.row==coord[1] & meta$spot.col==coord[2],])
    spot_info[spot_info$sample == sample & spot_info$spotid %in% spots, "cnv_cluster"] <- cnv_clusters[spot.id,]
  }
}   
#save spot data for all samples (overwrite the existing one)
write.table(spot_info, "", sep = "\t")

## CNV analysis

#metadata for samples
source("../utils/SamplesMetadata.R")

#load sub-spot data of all samples
spot_info <- read.table("")
cnv <- spot_info[!is.na(spot_info$cnv_cluster), ]

#Classify spots in TME, Cancer and Buffer depending on the estimate cluster
cnv$annot <- "TME"
cnv$annot[cnv$clusters == 3] <- "Buffer"
cnv$annot[cnv$clusters %in% c(1,2)] <- "Cancer"
#compute the percentage of TME, Buffer and Cancer spots in each of the CNV clusters
perc <- cnv %>% group_by(sample, cnv_cluster, annot) %>% 
  count() %>% group_by(sample, cnv_cluster) %>% mutate(percent=n/sum(n))  %>% 
  select(-n) %>% spread(annot, value = percent) %>% 
  mutate(BufferCancer = Buffer+Cancer) %>% filter(!is.na(cnv_cluster))
#select all CNV cluster with high content of Cancer (and Buffer) spots
perc_filter <- perc[(perc$Cancer > 0.8 | (perc$BufferCancer > 0.65) & perc$Cancer > 0.35),] %>% filter(!is.na(sample))
#select all spots and clusters that fulfill the above conditions
sub_cnv <- data.frame()
for (sample_s in unique(perc_filter$sample)) {
  for (cluster in filter(perc_filter, sample == sample_s)) {
    sub_cnv <- rbind(sub_cnv, cnv[cnv$sample==sample_s & cnv$cnv_cluster==cluster & cnv$clusters %in% c(1,2),])
  }
}
sub_cnv <- sub_cnv %>% filter(!is.na(sample))
colnames(sub_cnv) <- colnames(cnv)
sub_cnv$cnv_cluster <- factor(sub_cnv$cnv_cluster)
#compute the maximum activity difference (for each hallmark) across clones
test <- sub_cnv %>% group_by(sample, cnv_cluster) %>% summarise(across(H1:H13, mean)) %>% 
  group_by(sample) %>% summarise(minH1 = min(H1), maxH1=max(H1), diffH1=max(H1)-min(H1),
                                 minH2 = min(H2), maxH2=max(H2), diffH2=max(H2)-min(H2),
                                 minH3 = min(H3), maxH3=max(H3), diffH3=max(H3)-min(H3),
                                 minH4 = min(H4), maxH4=max(H4), diffH4=max(H4)-min(H4),
                                 minH5 = min(H5), maxH5=max(H5), diffH5=max(H5)-min(H5),
                                 minH6 = min(H6), maxH6=max(H6), diffH6=max(H6)-min(H6),
                                 minH7 = min(H7), maxH7=max(H7), diffH7=max(H7)-min(H7),
                                 minH8 = min(H8), maxH8=max(H8), diffH8=max(H8)-min(H8),
                                 minH9 = min(H9), maxH9=max(H9), diffH9=max(H9)-min(H9),
                                 minH10 = min(H10), maxH10=max(H10), diffH10=max(H10)-min(H10),
                                 minH11 = min(H11), maxH11=max(H11), diffH11=max(H11)-min(H11),
                                 minH12 = min(H12), maxH12=max(H12), diffH12=max(H12)-min(H12),
                                 minH13 = min(H13), maxH13=max(H13), diffH13=max(H13)-min(H13))
#create a data frame with these hallmark differences
df_diff <- data.frame(row.names = test$sample, test[, paste0("diffH", c(2,4,8,10,11,12))])
