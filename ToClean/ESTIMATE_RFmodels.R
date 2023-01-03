#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(caret)
df_all <- read.table("Desktop/df_all.txt", sep = "\t")

sample <- args[1]
hallmark <- args[2]

data <- df_all[df_all$sample == sample,]

fitControl <- trainControl(method = 'cv', number=10,summaryFunction=defaultSummary)

model = train(y=data[,hallmark], x=data[,c("estimate", "score")],
              method = "rf",
              preProc = c("center", "scale"),
              trControl = fitControl,
              metric='RMSE',
              trees=50,
              importance=T)


saveRDS(model, paste0("path/", sample, "_", hallmark, "_estimate_neigh.rds"))




