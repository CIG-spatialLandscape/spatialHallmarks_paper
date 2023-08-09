##################################################
## Project: Cancer Hallmarks
## Script purpose: Generate Random Forest model to predict a given neoplastic hallmark in a given sample
## Author: Sergi Cervilla* & Mustafa Sibai*
##################################################

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(caret)
library(ranger)
library(treeshap)
library(dplyr)

sample <- args[1] # sample name
target <- args[2] # target neoplastic Hallmark (one of H2, H4, H8, H9, H10, H11, H12)

radars <- read.table("", sep = "\t")

data <- select(radars, target, c(1:6))

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

