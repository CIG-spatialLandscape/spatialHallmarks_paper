library(ggplot2)

#### Load data ####
files <- list.files("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/output_df", full.names = F)
files <- stringr::str_remove(files, pattern = ".txt")
df_all <- lapply(files, function(x) {
  df <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/output_df/", x, ".txt"), sep = "\t")
  hallmark <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/objects_mts/", x, "_hallmarks.txt"), sep = "\t")
  hallmark <- hallmark[-1,]
  df[,paste0("H", 1:13)] <- hallmark
  return(df)
})
df_all <- do.call(rbind, df_all)

df_all$spotid <- sapply(rownames(df_all), function(spot){
  x <- strsplit(spot,  ".", fixed = T)
  paste0(x[[1]][1], ".", substr(x[[1]][2],1, 1))
})


# Scale estimate
df_all$estimate_scaled <- NA
for (sample in unique(df_all$sample)) {
  df_all$estimate_scaled[df_all$sample == sample] <- scale(df_all$estimate[df_all$sample == sample])
}

# Scale neighborhood
df_all$score_scaled <- NA
for (sample in unique(df_all$sample)) {
  df_all$score_scaled[df_all$sample == sample] <- scale(df_all$score[df_all$sample == sample])
}



################################################################################


###### lm: basal (whole) ######
lm_mat.basal <- data.frame(matrix(nrow = 0, ncol = 8))
for (sample in unique(df_all$sample)) {
  print(sample)
  sub_mat <-  df_all[df_all$sample == sample, ]
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
for (sample in unique(df_all$sample)) {
  print(sample)
  sub_mat <-  df_all[df_all$sample == sample, ]
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

## plots ##

# VolcanoPlot: estimate and significance
lm_mat.basal$meta <- paste0(lm_mat.basal$Sample, "~", lm_mat.basal$Hallmark)
lm_mat.basal$log <- -log10(lm_mat.basal$EstimatePvalue)
lm_mat.basal$log[lm_mat.basal$log == Inf] <- jitter(rep(320, length(lm_mat.basal$log[lm_mat.basal$log == Inf])))
lm_mat.basal$metric <-  abs(lm_mat.basal$Rsquared*lm_mat.basal$EstimateSlope)/lm_mat.basal$EstimateSlope

ggplot(lm_mat.basal, aes(x=metric, y=log)) +
  geom_point() + geom_text(data = lm_mat.basal[lm_mat.basal$log > 100,], aes(label=meta), check_overlap = T) +
  theme_classic() + labs(x="Rsquared*", y="-log10(p-value)") + theme(axis.text = element_text(size = 10),
                                                                     axis.title = element_text(size=15))
# VolcanoPlot: estimate and significance
lm_mat.1$meta <- paste0(lm_mat.1$Sample, "~", lm_mat.1$Hallmark)
lm_mat.1$log <- -log10(lm_mat.1$NeighborhoodPvalue)
lm_mat.1$log[lm_mat.1$log == Inf] <- jitter(rep(320, length(lm_mat.1$log[lm_mat.1$log == Inf])))
lm_mat.1$metric <-  abs(lm_mat.1$Rsquared*lm_mat.1$NeighborhoodSlope)/lm_mat.1$NeighborhoodSlope

ggplot(lm_mat.1, aes(x=NeighborhoodSlope, y=log)) +
  geom_point() + geom_text(data = lm_mat.1[lm_mat.1$log > 100,], aes(label=meta), check_overlap = T) +
  theme_classic() + labs(x="Rsquared*", y="-log10(p-value)") + theme(axis.text = element_text(size = 10),
                                                                     axis.title = element_text(size=15))
# Rsquared increase vs Rsquared basal
lm_mat.basal$diff <- lm_mat.1$Rsquared - lm_mat.basal$Rsquared
lm_mat.basal$metric <-  abs(lm_mat.basal$Rsquared*lm_mat.basal$EstimateSlope)/lm_mat.basal$EstimateSlope

ggplot(lm_mat.basal, aes(x=diff, y=Rsquared)) +
  geom_point() + geom_text(data = lm_mat.basal[lm_mat.basal$diff > 0.05,], aes(label=meta), check_overlap = T) +
  theme_classic() + labs(x="Rsquared increase", y="Rsquared (basal)") + theme(axis.text = element_text(size = 10),
                                                                              axis.title = element_text(size=15))



################################################################################

## Binary estimate


################################################################################
df_all$component <- NA
df_all$component[df_all$clusters %in% c(1,2)] <- "Cancer"
df_all$component[df_all$clusters %in% c(4,5)] <- "TME"
df_all$component <- factor(df_all$component, levels = c("Cancer", "TME"))
################################################################################


###### lm: basal (cluster) ######
lm_mat.basal.sub <- data.frame(matrix(nrow = 0, ncol = 8))
for (sample in unique(df_all$sample)) {
  print(sample)
  sub_mat <-  df_all[df_all$sample == sample & !is.na(df_all$component), ]
  for (hallmark in paste0("H", 1:13)) {
    f <- formula(paste0(hallmark, " ~  component"))
    lm_out <- lm(f, sub_mat)
    lm_summary <- summary(lm_out)
    lm_mat.basal.sub <- rbind(lm_mat.basal.sub, c(sample, hallmark, "all", lm_summary$adj.r.squared,
                                                  lm_summary$coefficients[2,1], lm_summary$coefficients[2,2], lm_summary$coefficients[2,3], lm_summary$coefficients[2,4]))
  }
}

colnames(lm_mat.basal.sub) <- c("Sample", "Hallmark", "Cluster", "Rsquared",
                                "EstimateSlope", "EstimateStdError", "EstimateTvalue", "EstimatePvalue")
lm_mat.basal.sub[,c(4:8)] <- apply(lm_mat.basal.sub[,c(4:8)], 2, as.numeric)

###### lm: estimate (binary) + neighborhood (1,2,4,5) ######
lm_mat.1.sub <- data.frame(matrix(nrow = 0, ncol = 12))
for (sample in unique(df_all$sample)) {
  print(sample)
  sub_mat <-  df_all[df_all$sample == sample, ]
  for (hallmark in paste0("H", 1:13)) {
    f <- formula(paste0(hallmark, " ~  component+score_scaled"))
    lm_out <- lm(f, sub_mat)
    lm_summary <- summary(lm_out)
    
    lm_mat.1.sub <- rbind(lm_mat.1.sub, c(sample, hallmark, "all", lm_summary$adj.r.squared,
                                          lm_summary$coefficients[2,1], lm_summary$coefficients[2,2], lm_summary$coefficients[2,3], lm_summary$coefficients[2,4],
                                          lm_summary$coefficients[3,1], lm_summary$coefficients[3,2], lm_summary$coefficients[3,3], lm_summary$coefficients[3,4]))
  }
}
colnames(lm_mat.1.sub) <- c("Sample", "Hallmark", "Cluster", "Rsquared",
                            "EstimateSlope", "EstimateStdError", "EstimateTvalue", "EstimatePvalue",
                            "NeighborhoodSlope", "NeighborhoodStdError", "NeighborhoodTvalue", "NeighborhoodPvalue")
lm_mat.1.sub[,c(4:12)] <- apply(lm_mat.1.sub[,c(4:12)], 2, as.numeric)

## plots ##




lm_mat.1.sub$meta <- paste0(lm_mat.1.sub$Sample, "~", lm_mat.1.sub$Hallmark)
lm_mat.1.sub$log <- -log10(lm_mat.1.sub$NeighborhoodPvalue)
lm_mat.1.sub$log[lm_mat.1.sub$log == Inf] <- jitter(rep(320, length(lm_mat.1.sub$log[lm_mat.1.sub$log == Inf])))
lm_mat.1.sub$metric <-  abs(lm_mat.1.sub$Rsquared*lm_mat.1.sub$NeighborhoodSlope)/lm_mat.1.sub$NeighborhoodSlope

ggplot(lm_mat.1.sub, aes(x=NeighborhoodSlope, y=log)) +
  geom_point() + geom_text(data = lm_mat.1.sub[lm_mat.1.sub$log > 100,], aes(label=meta), check_overlap = T) +
  theme_classic() + labs(x="Neighborhood slope", y="-log10(p-value)") + theme(axis.text = element_text(size = 10),
                                                                              axis.title = element_text(size=15))
## plots ##

# VolcanoPlot: estimate and significance
lm_mat.basal.sub$meta <- paste0(lm_mat.basal.sub$Sample, "~", lm_mat.basal.sub$Hallmark)
lm_mat.basal.sub$log <- -log10(lm_mat.basal.sub$EstimatePvalue)
lm_mat.basal.sub$log[lm_mat.basal.sub$log == Inf] <- jitter(rep(320, length(lm_mat.basal.sub$log[lm_mat.basal.sub$log == Inf])))
lm_mat.basal.sub$metric <-  abs(lm_mat.basal.sub$Rsquared*lm_mat.basal.sub$EstimateSlope)/lm_mat.basal.sub$EstimateSlope

ggplot(lm_mat.basal.sub, aes(x=metric, y=log)) +
  geom_point() + geom_text(data = lm_mat.basal.sub[lm_mat.basal.sub$log > 100,], aes(label=meta), check_overlap = T) +
  theme_classic() + labs(x="Rsquared*", y="-log10(p-value)") + theme(axis.text = element_text(size = 10),
                                                                     axis.title = element_text(size=15))


# VolcanoPlot: estimate and significance
lm_mat.1.sub$meta <- paste0(lm_mat.1.sub$Sample, "~", lm_mat.1.sub$Hallmark)
lm_mat.1.sub$log <- -log10(lm_mat.1.sub$NeighborhoodPvalue)
lm_mat.1.sub$log[lm_mat.1.sub$log == Inf] <- jitter(rep(320, length(lm_mat.1.sub$log[lm_mat.1.sub$log == Inf])))
lm_mat.1.sub$metric <-  abs(lm_mat.1.sub$Rsquared*lm_mat.1.sub$NeighborhoodSlope)/lm_mat.1.sub$NeighborhoodSlope

ggplot(lm_mat.1.sub, aes(x=NeighborhoodSlope, y=log)) +
  geom_point() + geom_text(data = lm_mat.1.sub[lm_mat.1.sub$log > 100,], aes(label=meta), check_overlap = T) +
  theme_classic() + labs(x="Rsquared*", y="-log10(p-value)") + theme(axis.text = element_text(size = 10),
                                                                     axis.title = element_text(size=15))
# Rsquared increase vs Rsquared basal
lm_mat.basal.sub$diff <- lm_mat.1.sub$Rsquared - lm_mat.basal.sub$Rsquared
lm_mat.basal.sub$metric <-  abs(lm_mat.basal.sub$Rsquared*lm_mat.basal.sub$EstimateSlope)/lm_mat.basal.sub$EstimateSlope

ggplot(lm_mat.basal.sub, aes(x=diff, y=Rsquared)) +
  geom_point() + geom_text(data = lm_mat.basal.sub[lm_mat.basal.sub$diff > 0.05,], aes(label=meta), check_overlap = T) +
  theme_classic() + labs(x="Rsquared increase", y="Rsquared (basal)") + theme(axis.text = element_text(size = 10),
                                                                              axis.title = element_text(size=15))

###### Scatter plots ######


ggplot(df_all[df_all$sample == "OVD1",], aes(x=estimate_scaled, y=H12)) +
  geom_point() + geom_smooth(method = "lm") + theme_classic() + ggtitle("OVD1 ~ H12") + theme(plot.title  = element_text(hjust = 0.5))
ggplot(df_all[df_all$sample == "C21",], aes(x=estimate_scaled, y=H5)) +
  geom_point() + geom_smooth(method = "lm") + theme_classic() + ggtitle("C21 ~ H5") + theme(plot.title  = element_text(hjust = 0.5))


ggplot(df_all[df_all$sample == "UKF242T",], aes(x=score_scaled, y=H5)) +
  geom_point() + geom_smooth(method = "lm") + theme_classic() + ggtitle("UKF242T ~ H5") + theme(plot.title  = element_text(hjust = 0.5))
ggplot(df_all[df_all$sample == "DU2",], aes(x=score_scaled, y=H10)) +
  geom_point() + geom_smooth(method = "lm") + theme_classic() + ggtitle("DU2 ~ H10") + theme(plot.title  = element_text(hjust = 0.5))


ggplot(df_all[df_all$sample == "UKF242T",], aes(x=estimate_scaled, y=H9)) +
  geom_point() + geom_smooth(method = "lm") + theme_classic() + ggtitle("UKF242T ~ H9") + theme(plot.title  = element_text(hjust = 0.5))
ggplot(df_all[df_all$sample == "UKF242T",], aes(x=score_scaled, y=H9)) +
  geom_point() + geom_smooth(method = "lm") + theme_classic() + ggtitle("UKF242T ~ H9") + theme(plot.title  = element_text(hjust = 0.5))

##
ggplot(df_all[df_all$sample == "HCC1T" & !is.na(df_all$component),], aes(x=score_scaled, y=H4)) +
  geom_point() + geom_smooth(method = "lm") + theme_classic() + ggtitle("HCC1T ~ H4") +
  theme(plot.title  = element_text(hjust = 0.5)) + facet_wrap(~component)
ggplot(df_all[df_all$sample == "C21" & !is.na(df_all$component),], aes(x=score_scaled, y=H4)) +
  geom_point() + geom_smooth(method = "lm") + theme_classic() + ggtitle("C21 ~ H4") +
  theme(plot.title  = element_text(hjust = 0.5)) + facet_wrap(~component)


ggplot(df_all[df_all$sample == "DU2" & !is.na(df_all$component),], aes(x=score_scaled, y=H6)) +
  geom_point() + geom_smooth(method = "lm") + theme_classic() + ggtitle("DU2 ~ H6") +
  theme(plot.title  = element_text(hjust = 0.5)) + facet_wrap(~component)
ggplot(df_all[df_all$sample == "C21" & !is.na(df_all$component),], aes(x=score_scaled, y=H11)) +
  geom_point() + geom_smooth(method = "lm") + theme_classic() + ggtitle("C21 ~ H11") +
  theme(plot.title  = element_text(hjust = 0.5)) + facet_wrap(~component)











## Whole


###### lm: basal (whole) ######
lm_mat.basal <- data.frame(matrix(nrow = 0, ncol = 8))
for (sample in unique(df_all$sample)) {
  print(sample)
  sub_mat <-  df_all[df_all$sample == sample, ]
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
for (sample in unique(df_all$sample)) {
  print(sample)
  sub_mat <-  df_all[df_all$sample == sample, ]
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
for (sample in unique(df_all$sample)) {
  print(sample)
  sub_mat <-  df_all[df_all$sample == sample & df_all$clusters %in% c(4,5), ]
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
for (sample in unique(df_all$sample)) {
  print(sample)
  sub_mat <-  df_all[df_all$sample == sample & df_all$clusters %in% c(4,5), ]
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
for (sample in unique(df_all$sample)) {
  print(sample)
  sub_mat <-  df_all[df_all$sample == sample & df_all$clusters %in% c(1,2), ]
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
for (sample in unique(df_all$sample)) {
  print(sample)
  sub_mat <-  df_all[df_all$sample == sample & df_all$clusters %in% c(1,2), ]
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
for (sample in unique(df_all$sample)) {
  print(sample)
  sub_mat <-  df_all[df_all$sample == sample & df_all$clusters %in% c(2,3,4), ]
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
for (sample in unique(df_all$sample)) {
  print(sample)
  sub_mat <-  df_all[df_all$sample == sample & df_all$clusters %in% c(2,3,4), ]
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


# 
# ############################
# ###### lm: basal (whole) ######
# lm_mat.basal <- data.frame(matrix(nrow = 0, ncol = 8))
# for (sample in unique(df_all$sample)) {
#   print(sample)
#   sub_mat <-  df_all[df_all$sample == sample, ]
#   for (hallmark in paste0("H", 1:13)) {
#     f <- formula(paste0(hallmark, " ~  estimate"))
#     lm_out <- lm(f, sub_mat)
#     lm_summary <- summary(lm_out)
#     
#     lm_mat.basal <- rbind(lm_mat.basal, c(sample, hallmark, "all", lm_summary$adj.r.squared,
#                                           lm_summary$coefficients[2,1], lm_summary$coefficients[2,2], lm_summary$coefficients[2,3], lm_summary$coefficients[2,4]))
#   }
# }
# 
# colnames(lm_mat.basal) <- c("Sample", "Hallmark", "Cluster", "Rsquared",
#                             "EstimateSlope", "EstimateStdError", "EstimateTvalue", "EstimatePvalue")
# 
# lm_mat.basal[,c(4:8)] <- apply(lm_mat.basal[,c(4:8)], 2, as.numeric)


# ###### lm: residuals (whole) ######
# lm_mat.residuals <- data.frame(matrix(nrow = 0, ncol = 8))
# for (sample in unique(df_all$sample)) {
#   print(sample)
#   sub_mat <-  df_all[df_all$sample == sample, ]
#   #compute residuals
#   residuals_out <- lm(estimate~score, sub_mat)
#   sub_mat$residuals <- residuals_out$residuals
#   for (hallmark in paste0("H", 1:13)) {
#     
#     f <- formula(paste0(hallmark, " ~  estimate + residuals"))
#     lm_out <- lm(f, sub_mat)
#     lm_summary <- summary(lm_out)
#     
#     lm_mat.residuals <- rbind(lm_mat.residuals, c(sample, hallmark, "residuals", lm_summary$adj.r.squared,
#                                           lm_summary$coefficients[2,1], lm_summary$coefficients[2,2], lm_summary$coefficients[2,3], lm_summary$coefficients[2,4]))
#   }
# }
# 
# colnames(lm_mat.residuals) <- c("Sample", "Hallmark", "Cluster", "Rsquared",
#                             "EstimateSlope", "EstimateStdError", "EstimateTvalue", "EstimatePvalue")
# 
# lm_mat.residuals[,c(4:8)] <- apply(lm_mat.residuals[,c(4:8)], 2, as.numeric)

# ###### lm: estimate + neighborhood (whole) ######
# lm_mat.1 <- data.frame(matrix(nrow = 0, ncol = 12))
# for (sample in unique(df_all$sample)) {
#   print(sample)
#   sub_mat <-  df_all[df_all$sample == sample, ]
#   for (hallmark in paste0("H", 1:13)) {
#     f <- formula(paste0(hallmark, " ~  estimate_scaled+score_scaled"))
#     lm_out <- lm(f, sub_mat)
#     lm_summary <- summary(lm_out)
#     
#     lm_mat.1 <- rbind(lm_mat.1, c(sample, hallmark, "score", lm_summary$adj.r.squared,
#                                   lm_summary$coefficients[2,1], lm_summary$coefficients[2,2], lm_summary$coefficients[2,3], lm_summary$coefficients[2,4],
#                                   lm_summary$coefficients[3,1], lm_summary$coefficients[3,2], lm_summary$coefficients[3,3], lm_summary$coefficients[3,4]))
#   }
# }
# colnames(lm_mat.1) <- c("Sample", "Hallmark", "Cluster", "Rsquared",
#                         "EstimateSlope", "EstimateStdError", "EstimateTvalue", "EstimatePvalue",
#                         "NeighborhoodSlope", "NeighborhoodStdError", "NeighborhoodTvalue", "NeighborhoodPvalue")
# 
# 
# lm_mat.1[,c(4:12)] <- apply(lm_mat.1[,c(4:12)], 2, as.numeric)
# 
# df <- rbind(lm_mat.basal, lm_mat.residuals, lm_mat.1[,1:8])
# ggplot(df, aes(x=Cluster, y=Rsquared)) + geom_boxplot()
