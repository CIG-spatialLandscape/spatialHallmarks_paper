##################################################
## Project: Cancer Hallmarks
## Script purpose: Generate Random Forest model to predict a given TME hallmark in a given sample
## Author: Sergi Cervilla* & Mustafa Sibai*
##################################################

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(caret)
library(ranger)
library(treeshap)
library(dplyr)

sample <- args[1] # sample name
target <- args[2] # target TME hallmark (one of H1, H3, H5, H6, H7, H13)

radars <- read.table("", sep = "\t")

data <- dplyr::select(radars, target, c(1:7))

inTrain <- createDataPartition(y = data[,target], p = 0.8, list = FALSE)

training <- data[inTrain,]

testing <- data[-inTrain,]


set.seed(123)

f <- formula(paste0(target, "~."))
model <- ranger::ranger(f, data = training, importance = "impurity", num.trees = 500, scale.permutation.importance = T)


model_unified <- ranger.unify(model, data)


treeshap_res <- treeshap(model_unified, testing)

saveRDS(model)
saveRDS(treeshap_res)
