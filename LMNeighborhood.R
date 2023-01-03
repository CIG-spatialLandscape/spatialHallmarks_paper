##################################################
## Project: Cancer Hallmarks
## Script purpose: Linear models to predict hallmark activity through estimate and neighborhood scores
## Date: 22/12/2022
## Author: Sergi Cervilla & Mustafa Sibai
##################################################

#### Load data ####
#extract all samples names from files
files <- list.files("", full.names = F)
files <- stringr::str_remove(files, pattern = ".txt")
#unify files and samples
spot_info <- lapply(files, function(x) {
  #spot data and neighborhoodscore
  df <- read.table("", sep = "\t")
  #hallmark activity for each sub-spot
  hallmark <- read.table("", sep = "\t")
  hallmark <- hallmark[-1,]
  df[,paste0("H", 1:13)] <- hallmark
  return(df)
})
spot_info <- do.call(rbind, spot_info)

spot_info$spotid <- sapply(rownames(spot_info), function(spot){
  x <- strsplit(spot,  ".", fixed = T)
  paste0(x[[1]][1], ".", substr(x[[1]][2],1, 1))
})


# Scale estimate
spot_info$estimate_scaled <- NA
for (sample in unique(spot_info$sample)) {
  spot_info$estimate_scaled[spot_info$sample == sample] <- scale(spot_info$estimate[spot_info$sample == sample])
}

# Scale neighborhood
spot_info$score_scaled <- NA
for (sample in unique(spot_info$sample)) {
  spot_info$score_scaled[spot_info$sample == sample] <- scale(spot_info$score[spot_info$sample == sample])
}


# Compute linear models

## Whole
###### lm: basal (whole) ######
lm_mat.basal <- data.frame(matrix(nrow = 0, ncol = 8))
for (sample in unique(spot_info$sample)) {
  print(sample)
  sub_mat <-  spot_info[spot_info$sample == sample, ]
  for (hallmark in paste0("H", 1:13)) {
    f <- formula(paste0(hallmark, " ~  estimate"))
    lm_out <- lm(f, sub_mat)
    lm_summary <- summary(lm_out)
    
    lm_mat.basal <- rbind(lm_mat.basal, c(sample, hallmark, "all", lm_summary$adj.r.squared,
                                          lm_summary$coefficients[2,1], lm_summary$coefficients[2,2], lm_summary$coefficients[2,3], lm_summary$coefficients[2,4]))
  }
}

colnames(lm_mat.basal) <- c("Sample", "Hallmark", "Cluster", "Rsquared",
                            "EstimateSlope", "EstimateStdError", "EstimateTvalue", "EstimatePvalue")

lm_mat.basal[,c(4:8)] <- apply(lm_mat.basal[,c(4:8)], 2, as.numeric)

###### lm: estimate + neighborhood (whole) ######
lm_mat.1 <- data.frame(matrix(nrow = 0, ncol = 12))
for (sample in unique(spot_info$sample)) {
  print(sample)
  sub_mat <-  spot_info[spot_info$sample == sample, ]
  for (hallmark in paste0("H", 1:13)) {
    f <- formula(paste0(hallmark, " ~  estimate_scaled+score_scaled"))
    lm_out <- lm(f, sub_mat)
    lm_summary <- summary(lm_out)
    
    lm_mat.1 <- rbind(lm_mat.1, c(sample, hallmark, "all", lm_summary$adj.r.squared,
                                  lm_summary$coefficients[2,1], lm_summary$coefficients[2,2], lm_summary$coefficients[2,3], lm_summary$coefficients[2,4],
                                  lm_summary$coefficients[3,1], lm_summary$coefficients[3,2], lm_summary$coefficients[3,3], lm_summary$coefficients[3,4]))
  }
}
colnames(lm_mat.1) <- c("Sample", "Hallmark", "Cluster", "Rsquared",
                        "EstimateSlope", "EstimateStdError", "EstimateTvalue", "EstimatePvalue",
                        "NeighborhoodSlope", "NeighborhoodStdError", "NeighborhoodTvalue", "NeighborhoodPvalue")


lm_mat.1[,c(4:12)] <- apply(lm_mat.1[,c(4:12)], 2, as.numeric)


## TME


###### lm: basal (whole) ######
lm_mat.basal.TME <- data.frame(matrix(nrow = 0, ncol = 8))
for (sample in unique(spot_info$sample)) {
  print(sample)
  sub_mat <-  spot_info[spot_info$sample == sample & spot_info$clusters %in% c(4,5), ]
  for (hallmark in paste0("H", 1:13)) {
    f <- formula(paste0(hallmark, " ~  estimate"))
    lm_out <- lm(f, sub_mat)
    lm_summary <- summary(lm_out)
    
    lm_mat.basal.TME <- rbind(lm_mat.basal.TME, c(sample, hallmark, "all", lm_summary$adj.r.squared,
                                                  lm_summary$coefficients[2,1], lm_summary$coefficients[2,2], lm_summary$coefficients[2,3], lm_summary$coefficients[2,4]))
  }
}

colnames(lm_mat.basal.TME) <- c("Sample", "Hallmark", "Cluster", "Rsquared",
                                "EstimateSlope", "EstimateStdError", "EstimateTvalue", "EstimatePvalue")

lm_mat.basal.TME[,c(4:8)] <- apply(lm_mat.basal.TME[,c(4:8)], 2, as.numeric)

###### lm: estimate + neighborhood (whole) ######
lm_mat.1.TME <- data.frame(matrix(nrow = 0, ncol = 12))
for (sample in unique(spot_info$sample)) {
  print(sample)
  sub_mat <-  spot_info[spot_info$sample == sample & spot_info$clusters %in% c(4,5), ]
  for (hallmark in paste0("H", 1:13)) {
    f <- formula(paste0(hallmark, " ~  estimate_scaled+score_scaled"))
    lm_out <- lm(f, sub_mat)
    lm_summary <- summary(lm_out)
    
    lm_mat.1.TME <- rbind(lm_mat.1.TME, c(sample, hallmark, "all", lm_summary$adj.r.squared,
                                          lm_summary$coefficients[2,1], lm_summary$coefficients[2,2], lm_summary$coefficients[2,3], lm_summary$coefficients[2,4],
                                          lm_summary$coefficients[3,1], lm_summary$coefficients[3,2], lm_summary$coefficients[3,3], lm_summary$coefficients[3,4]))
  }
}
colnames(lm_mat.1.TME) <- c("Sample", "Hallmark", "Cluster", "Rsquared",
                            "EstimateSlope", "EstimateStdError", "EstimateTvalue", "EstimatePvalue",
                            "NeighborhoodSlope", "NeighborhoodStdError", "NeighborhoodTvalue", "NeighborhoodPvalue")


lm_mat.1.TME[,c(4:12)] <- apply(lm_mat.1.TME[,c(4:12)], 2, as.numeric)




## Cancer

###### lm: basal (whole) ######
lm_mat.basal.Cancer <- data.frame(matrix(nrow = 0, ncol = 8))
for (sample in unique(spot_info$sample)) {
  print(sample)
  sub_mat <-  spot_info[spot_info$sample == sample & spot_info$clusters %in% c(1,2), ]
  for (hallmark in paste0("H", 1:13)) {
    f <- formula(paste0(hallmark, " ~  estimate"))
    lm_out <- lm(f, sub_mat)
    lm_summary <- summary(lm_out)
    
    lm_mat.basal.Cancer <- rbind(lm_mat.basal.Cancer, c(sample, hallmark, "all", lm_summary$adj.r.squared,
                                                        lm_summary$coefficients[2,1], lm_summary$coefficients[2,2], lm_summary$coefficients[2,3], lm_summary$coefficients[2,4]))
  }
}

colnames(lm_mat.basal.Cancer) <- c("Sample", "Hallmark", "Cluster", "Rsquared",
                                   "EstimateSlope", "EstimateStdError", "EstimateTvalue", "EstimatePvalue")

lm_mat.basal.Cancer[,c(4:8)] <- apply(lm_mat.basal.Cancer[,c(4:8)], 2, as.numeric)

###### lm: estimate + neighborhood (whole) ######
lm_mat.1.Cancer <- data.frame(matrix(nrow = 0, ncol = 12))
for (sample in unique(spot_info$sample)) {
  print(sample)
  sub_mat <-  spot_info[spot_info$sample == sample & spot_info$clusters %in% c(1,2), ]
  for (hallmark in paste0("H", 1:13)) {
    f <- formula(paste0(hallmark, " ~  estimate_scaled+score_scaled"))
    lm_out <- lm(f, sub_mat)
    lm_summary <- summary(lm_out)
    
    lm_mat.1.Cancer <- rbind(lm_mat.1.Cancer, c(sample, hallmark, "all", lm_summary$adj.r.squared,
                                                lm_summary$coefficients[2,1], lm_summary$coefficients[2,2], lm_summary$coefficients[2,3], lm_summary$coefficients[2,4],
                                                lm_summary$coefficients[3,1], lm_summary$coefficients[3,2], lm_summary$coefficients[3,3], lm_summary$coefficients[3,4]))
  }
}
colnames(lm_mat.1.Cancer) <- c("Sample", "Hallmark", "Cluster", "Rsquared",
                               "EstimateSlope", "EstimateStdError", "EstimateTvalue", "EstimatePvalue",
                               "NeighborhoodSlope", "NeighborhoodStdError", "NeighborhoodTvalue", "NeighborhoodPvalue")


lm_mat.1.Cancer[,c(4:12)] <- apply(lm_mat.1.Cancer[,c(4:12)], 2, as.numeric)




## Buffer

###### lm: basal (whole) ######
lm_mat.basal.Buffer <- data.frame(matrix(nrow = 0, ncol = 8))
for (sample in unique(spot_info$sample)) {
  print(sample)
  sub_mat <-  spot_info[spot_info$sample == sample & spot_info$clusters %in% c(2,3,4), ]
  for (hallmark in paste0("H", 1:13)) {
    f <- formula(paste0(hallmark, " ~  estimate"))
    lm_out <- lm(f, sub_mat)
    lm_summary <- summary(lm_out)
    
    lm_mat.basal.Buffer <- rbind(lm_mat.basal.Buffer, c(sample, hallmark, "all", lm_summary$adj.r.squared,
                                                        lm_summary$coefficients[2,1], lm_summary$coefficients[2,2], lm_summary$coefficients[2,3], lm_summary$coefficients[2,4]))
  }
}

colnames(lm_mat.basal.Buffer) <- c("Sample", "Hallmark", "Cluster", "Rsquared",
                                   "EstimateSlope", "EstimateStdError", "EstimateTvalue", "EstimatePvalue")

lm_mat.basal.Buffer[,c(4:8)] <- apply(lm_mat.basal.Buffer[,c(4:8)], 2, as.numeric)

###### lm: estimate + neighborhood (whole) ######
lm_mat.1.Buffer <- data.frame(matrix(nrow = 0, ncol = 12))
for (sample in unique(spot_info$sample)) {
  print(sample)
  sub_mat <-  spot_info[spot_info$sample == sample & spot_info$clusters %in% c(2,3,4), ]
  for (hallmark in paste0("H", 1:13)) {
    f <- formula(paste0(hallmark, " ~  estimate_scaled+score_scaled"))
    lm_out <- lm(f, sub_mat)
    lm_summary <- summary(lm_out)
    
    lm_mat.1.Buffer <- rbind(lm_mat.1.Buffer, c(sample, hallmark, "all", lm_summary$adj.r.squared,
                                                lm_summary$coefficients[2,1], lm_summary$coefficients[2,2], lm_summary$coefficients[2,3], lm_summary$coefficients[2,4],
                                                lm_summary$coefficients[3,1], lm_summary$coefficients[3,2], lm_summary$coefficients[3,3], lm_summary$coefficients[3,4]))
  }
}
colnames(lm_mat.1.Buffer) <- c("Sample", "Hallmark", "Cluster", "Rsquared",
                               "EstimateSlope", "EstimateStdError", "EstimateTvalue", "EstimatePvalue",
                               "NeighborhoodSlope", "NeighborhoodStdError", "NeighborhoodTvalue", "NeighborhoodPvalue")


lm_mat.1.Buffer[,c(4:12)] <- apply(lm_mat.1.Buffer[,c(4:12)], 2, as.numeric)



#merge in one table all important information
lm_mat.basal$diff <- lm_mat.1$Rsquared - lm_mat.basal$Rsquared
lm_mat.basal$RsquaredCancer <- lm_mat.basal.Cancer$Rsquared
lm_mat.basal$diffCancer <- lm_mat.1.Cancer$Rsquared - lm_mat.basal.Cancer$Rsquared
lm_mat.basal$RsquaredTME <- lm_mat.basal.TME$Rsquared
lm_mat.basal$diffTME <- lm_mat.1.TME$Rsquared - lm_mat.basal.TME$Rsquared
lm_mat.basal$RsquaredBuffer <- lm_mat.basal.Buffer$Rsquared
lm_mat.basal$diffBuffer <- lm_mat.1.Buffer$Rsquared - lm_mat.basal.Buffer$Rsquared


# Define CancerType groups

lm_mat.basal$CancerType <- NA
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("Acinar", "IC", "PC1", "PC2")] <- "Prostate"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("Breast", "BreastA", "Ductal", "DuctalFFPE", "TNBCA", "M1", "M2", "M3", "M4")] <- "Breast"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("C20", "C21", "C34", "C51", "C7")] <- "Kidney"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("cHC1T", "HCC1T", "HCC2T", "HCC5D", "ICC1L")] <- "Liver"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("Colorectal", "CRC1", "CRC2", "Intestine", "Co1", "Co2", "Co3", "Co4")] <- "Colorectal"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("DU2", "DU3", "DU12", "DU8", "DU13")] <- "Bladder"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("OV4A", "OVD1", "Ovarian", "CUP295", "OVFFPE")] <- "Ovarian"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("P259_H2A2", "P264", "P270", "P288", "P306")] <- "Pancreas"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("Glioblastoma", "UKF242T", "UKF260T", "UKF269T", "UKF275T")] <- "Glioblastoma"
lm_mat.basal$CancerType[lm_mat.basal$Sample %in% c("P1","P3", "P4", "P5", "P6", "P8")] <- "Lung"


