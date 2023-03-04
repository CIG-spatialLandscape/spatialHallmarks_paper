##################################################
## Project: Cancer Hallmarks
## Script purpose: Compute the pathway scores in an enhanced ST object
## Date: 22/12/2022
## Author: Sergi Cervilla * & Mustafa Sibai *
##################################################

#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = T)

library(Seurat)
library(STutility)
library(dplyr)
library(RColorBrewer)
library(NNLM)
library(viridis)
library(SingleCellExperiment)
library(BayesSpace)

sample <- args
#load enhanced Seurat object
STobject_enhanced <- readRDS("")
#set SCT assay as default
DefaultAssay(STobject_enhanced) <- "SCT"

#Concatenate all pathway names of each hallmark
path.names <- c()
for (hallmark in paste0("H", 1:13)) {
  #load list of pathways for a given hallmark
  H_paths <- readRDS("")
  #pathways that at least have 2 genes in the Seurat object
  H_paths.2 <- H_paths
  for (path in 1:length(H_paths)) {
    #Check how many genes of that pathway are present in the Seurat object
    if(sum(H_paths[[path]] %in% rownames(STobject_enhanced)) < 1) {
      print(names(H_paths[path]))
      #remove pathways with less than 2 genes
      H_paths.2[names(H_paths[path])] <- NULL
    }
  }
  #vector of names with pathway names with more than 1 gene (from a hallmark)
  path.names <- c(path.names, paste0(hallmark, "_", names(H_paths.2)))
  
  #compute Pathway activity through AddModuleScore
  STobject_enhanced <- AddModuleScore(
    STobject_enhanced,
    H_paths.2,
    pool = NULL,
    nbin = 24,
    ctrl = 100,
    k = FALSE,
    assay = NULL,
    name = paste0(hallmark, "_"),
    seed = 200,
    search = F
  )
}


#obtain estimate cluster and hallmark activities
metadata <- STobject_enhanced@meta.data[,c("estimate.cluster", grep("H", colnames(STobject_enhanced@meta.data),value = T))]
#replace pathway names to their corresponding name
colnames(metadata)[15:ncol(metadata)] <- path.names
#save the matrix
write.table(metadata, "", sep= "\t", quote = F)
