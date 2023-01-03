#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(caret)
library(randomForest)
df_tme <- read.table("Desktop/df_tme.txt", sep = "\t")

sample <- args[1]
hallmark <- args[2]

data <- df_tme[df_tme$sample == sample,]

fitControl <- trainControl(method = 'cv', number=10,summaryFunction=defaultSummary)

model = train(y=data[,hallmark], x=data[,c(19:24)],
              method = "rf",
              preProc = c("center", "scale"),
              trControl = fitControl,
              metric='RMSE',
              trees=ntrees,
              importance=T)


saveRDS(model, paste0("path/", sample, "_", hallmark, ".rds"))




