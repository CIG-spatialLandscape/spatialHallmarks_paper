##################################################
## Project: Cancer Hallmarks
## Script purpose: Impute genes at sub-spot resolution
## Date: 22/12/2022
## Author: Sergi Cervilla & Mustafa Sibai
##################################################

#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = T)

library(BayesSpace)

sample <- args
#load single cell experiment object at spot level (with gene expression)
sce <- readRDS("")
#load single cell experiment object at sub-spot level (without gene expression)
sce_enhanced <- readRDS("")

#Impute genes at sub-spot resolution
sce_enhanced <- enhanceFeatures(sce_enhanced, sce,
                                model="xgboost",
                                nrounds=100,
                                assay.type = "SCT")
#save file
saveRDS(sce_enhanced, "")
