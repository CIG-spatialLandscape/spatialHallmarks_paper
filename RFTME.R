##################################################
## Project: Cancer Hallmarks
## Script purpose: Generate Random Forest model to predict a given TME hallmark in a given sample
## Date: 22/12/2022
## Author: Sergi Cervilla * & Mustafa Sibai *
##################################################

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(caret)
library(randomForest)
library(ranger)
library(treeshap)
library(dplyr)

#data frame containing hallmark activity and radar scores
df_tme <- read.table("")

sample <- args[1] #first argument will be sample name
hallmark <- args[2] #second argument will be hallmark to predict (e.g. H1)

#subset by sample
data <- df_tme[df_tme$sample == sample,]
#select predictor (a TME hallmark) and Radar scores for Cancer hallmarks (predictors)
data <- select(data, hallmark, c(18:23))

#split into training a test data
inTrain <- createDataPartition(y = data[,hallmark], p = 0.8, list = FALSE)
training <- data[inTrain,]
testing <- data[-inTrain,]


set.seed(123)
#formula for the Random Forest Regressor
f <- formula(paste0(hallmark, "~."))
#model for the Random Forest Regressor
model <- ranger::ranger(f, data = training, importance = "impurity", num.trees = 500, scale.permutation.importance = T)
model_unified <- ranger.unify(model, data)
#compute shapley additive values (importance and dependency of each predictor)
treeshap_res <- treeshap(model_unified, testing)
#save model
saveRDS(model, "")
#save shapley results
saveRDS(treeshap_res, "")

