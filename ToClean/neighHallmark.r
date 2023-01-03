new_labels <- function(original, order){
  labels <- original
  #Cancer (lowest)
  labels[original==order[1]] <- "1"
  labels[original==order[2]] <- "2"
  labels[original==order[3]] <- "3"
  labels[original==order[4]] <- "4"
  #TME (highest)
  labels[original==order[5]] <- "5"
  return(labels)
}


select_neighbours <- function(spot, threshold, distances) {
  return(grep(x = colnames(distances)[distances[spot, ] < threshold], pattern = spot, value = T, invert = T))
}



STobject <- readRDS("Desktop/enhanced/Seurat/Colorectal_enhanced_int.rds")
STobject <- STobject[,-1]
coord <- STobject@images[[1]]@coordinates
hallmarks <- STobject@meta.data[, paste0("H", 1:13)]
sample <- "Colorectal"
rm(STobject)
gc()
estimate <- read.table(paste0("Desktop/estimate/output/output_", sample, ".gct"), sep = "\t", skip = 2, header = T)
estimate <- t(estimate[3, 3:ncol(estimate)])

##### fix coordinates
for (subspot in rownames(coord)) {
  realrow <-  coord[paste0(substr(subspot, 0, nchar(subspot)-2), ".5"), "row"]
  realrow <- (realrow+1)*100
  realcol <- trunc(coord[paste0(substr(subspot, 0, nchar(subspot)-2), ".3"), "col"])
  realcol <- (realcol+1)*100
  n <- substr(subspot, nchar(subspot), nchar(subspot))
  
  factor <- 55/3
  if (n == 1) {
    realrow <- realrow + factor
    realcol <- realcol + factor
  }  else if (n == 2) {
    realrow <- realrow + factor
    realcol <- realcol - factor
  }  else if (n == 3) {
    realrow <- realrow - factor
    realcol <- realcol + factor
  }  else if (n == 4) {
    realrow <- realrow - factor
    realcol <- realcol - factor
  }  else if (n == 5) {
    realrow <- realrow
    realcol <- realcol + 2 * factor
  }  else if (n == 6) {
    realrow <- realrow
    realcol <- realcol - 2 * factor
  }
  coord[subspot, c("realrow", "realcol")] <- c(realrow, realcol)
}

##### select neigh
distances <- as.matrix(dist(coord[,c("realrow", "realcol")]))


score <- sapply(rownames(coord), function(spot) {
  neigh <- select_neighbours(spot, 200, distances)
  score <- mean(estimate[neigh,]/distances[spot, neigh])
})



#estimate <- scale(estimate)
df <- data.frame(score = score)
df <- cbind(df, scale(hallmarks))
labels <- kmeans(estimate, centers = 5)$cluster

labels_order <- order(sapply(1:5, function(x){
  mean(estimate[labels==x])
}))

new_labels <- function(original, order){
  labels <- original
  #Cancer (lowest)
  labels[original==order[1]] <- "1"
  labels[original==order[2]] <- "2"
  labels[original==order[3]] <- "3"
  labels[original==order[4]] <- "4"
  #TME (highest)
  labels[original==order[5]] <- "5"
  return(labels)
}

df$clusters <- new_labels(labels,labels_order)
df$estimate <- estimate


ggplot(df, aes(x=score, y=H12)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~clusters)

  ggplot(df, aes(x=score, y=H7)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~clusters)


cor(df$score[df$clusters==1], df$estimate[df$clusters==1])  
cor(df$score[df$clusters==2], df$estimate[df$clusters==2])  
cor(df$score[df$clusters==3], df$estimate[df$clusters==3])  
cor(df$score[df$clusters==4], df$estimate[df$clusters==4])  
cor(df$score[df$clusters==5], df$estimate[df$clusters==5])  
ggplot(df, aes(x=score, y=estimate)) + geom_point() + facet_wrap(~clusters)


########
  sce <- readRDS("Desktop/enhanced/enhanced_sce/Colorectal_sce_enhanced.rds") 
  
  source("Desktop/IJC/datasets/IGTP/figuresPaper/scripts/utilities/BayesSpace_functions.r")
  sce <- sce[,-1]
  sce$score <- score*estimate**2
  v <- .make_triangle_subspots(colData(sce), fill = "score")
  p <- ggplot() +  geom_polygon(data=v,
                                aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~(X3))) + theme_void()
  p
  sce$estimate <- estimate
  v <- .make_triangle_subspots(colData(sce), fill = "estimate")
  p2 <- ggplot() +  geom_polygon(data=v,
                                aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~(X3))) 
p+p2  

write.table(df, "Desktop/OV4A_neigh.txt", sep = "\t")



df1 <- read.table( "Desktop/OV4A_neigh.txt", sep = "\t")
df$sample <- "Colorectal"
df1$sample <- "OV4A"


df_all <- rbind(df, df1)
ggplot(df_all, aes(x=score, y=H12, col = sample)) + geom_point()



# Deal with batch effects

## Experiment 1 - scale score (s)

df_scaled <- df
df_scaled$score <- scale(df$score)

df1_scaled <- df1
df1_scaled$score <- scale(df1$score)


df_all_scaled <- rbind(df_scaled, df1_scaled)
ggplot(df_all_scaled, aes(x=score, y=H12, col = sample)) + geom_point()


## Experiment 1 - scale score (0-1)

df_scaled <- df
df_scaled$score <- scale(df$score)

df1_scaled <- df1
df1_scaled$score <- scale(df1$score)


df_all_scaled <- rbind(df_scaled, df1_scaled)
ggplot(df_all_scaled, aes(x=score, y=H12, col = sample)) + geom_point()


################################################################################

samples <- unique(sapply(list.files("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/objects_mts/"), function(x){
  strsplit(x, split = "_")[[1]][1]
}))
samples[32] <- "P259_H2A2"

samples <- c("P3","P4", "P6", "PC1", "PC2")
samples <- c("Co1", "Co2", "Co3", "Co4", "M1", "M2", "M3", "M4", "P1", "P5", "P8")
samples <- "P7"
for (sample in samples){
  print(sample)
  hallmarks <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/objects_mts/", sample, "_hallmarks.txt"), sep = "\t", header = T)
  hallmarks <- hallmarks[rownames(hallmarks)!="subspot_1.1",]
  coord <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/objects_mts/", sample, "_coords.txt"), sep = "\t", header = T)
  coord <- coord[rownames(coord)!="subspot_1.1",]
  estimate <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/estimate/output/output_", sample, ".gct"), sep = "\t", skip = 2, header = T)
  estimate <- t(estimate[3, 3:ncol(estimate)])
  
  ##### fix coordinates
  for (subspot in rownames(coord)) {
    realrow <-  coord[paste0(substr(subspot, 0, nchar(subspot)-2), ".5"), "row"]
    realrow <- (realrow+1)*100
    realcol <- trunc(coord[paste0(substr(subspot, 0, nchar(subspot)-2), ".3"), "col"])
    realcol <- (realcol+1)*100
    n <- substr(subspot, nchar(subspot), nchar(subspot))
    
    factor <- 55/3
    if (n == 1) {
      realrow <- realrow + factor
      realcol <- realcol + factor
    }  else if (n == 2) {
      realrow <- realrow + factor
      realcol <- realcol - factor
    }  else if (n == 3) {
      realrow <- realrow - factor
      realcol <- realcol + factor
    }  else if (n == 4) {
      realrow <- realrow - factor
      realcol <- realcol - factor
    }  else if (n == 5) {
      realrow <- realrow
      realcol <- realcol + 2 * factor
    }  else if (n == 6) {
      realrow <- realrow
      realcol <- realcol - 2 * factor
    }
    coord[subspot, c("realrow", "realcol")] <- c(realrow, realcol)
  }
  
  ##### select neigh
  distances <- as.matrix(dist(coord[,c("realrow", "realcol")]))
  gc()
  
  score <- sapply(rownames(hallmarks), function(spot) {
    neigh <- select_neighbours(spot, 200, distances)
    score <- mean(estimate[neigh,]/distances[spot, neigh])
  })
  
  df <- data.frame(score = score)
  df <- cbind(df, scale(hallmarks))
  labels <- kmeans(estimate, centers = 5)$cluster
  
  labels_order <- order(sapply(1:5, function(x){
    mean(estimate[labels==x])
  }))
  
  df$clusters <- new_labels(labels,labels_order)
  df$estimate <- estimate
  df$sample <- sample
  write.table(df, paste0("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/output_df/", sample, ".txt"), sep = "\t")
  rm(estimate, df, hallmarks, distances, coord)
  gc()
}





files <- list.files("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/output_df", full.names = T)


df_all <- lapply(files, function(x) {
  read.table(x, sep = "\t")
})
df_all <- do.call(rbind, df_all)

ggplot(df_all, aes(x=score, y=H12, col=sample)) + geom_point()



#####3

ggplot(df_all, aes(x=score, y=sample)) + geom_density_ridges()

ggplot(df_all, aes(x=estimate, y=sample)) + geom_density_ridges()


df_all$estimate_scaled <- NA
for (sample in unique(df_all$sample)) {
  df_all$estimate_scaled[df_all$sample == sample] <- scale(df_all$estimate[df_all$sample == sample])
}
ggplot(df_all, aes(x=estimate_scaled, y=sample)) + geom_density_ridges()

df_all$score_scaled <- NA
for (sample in unique(df_all$sample)) {
  df_all$score_scaled[df_all$sample == sample] <- scale(df_all$score[df_all$sample == sample])
}
ggplot(df_all, aes(x=score_scaled, y=sample)) + geom_density_ridges()

ggplot(df_all, aes(x=score_scaled, y=H2, col=sample)) + geom_point()
ggplot(df_all, aes(x=score, y=H12, col=sample)) + geom_smooth(method = "lm")

ggplot(df_all, aes(x=score, y=H2)) + geom_smooth(method = "lm") + facet_wrap(~sample, scales = "free")




ggplot(df_all[df_all$clusters==1,], aes(x=score, y=H2)) + geom_smooth(method = "lm") + facet_wrap(~sample, scales = "free")




ggplot(df_all[df_all$sample=="CRC1",], aes(x=score, y=H12, col=estimate)) + geom_point()  + geom_smooth(method = "lm")

ggplot(df_all[abs(df_all$score_scaled) < 5, ], aes(x=score_scaled, y=H12)) + geom_smooth(method = "lm") + facet_wrap(~clusters)
ggplot(df_all[abs(df_all$score_scaled) < 5, ], aes(x=score_scaled, y=H12)) + geom_smooth(method = "lm")

library(lme4)
model <- lmer(H12 ~ score_scaled + (1 + score_scaled | sample), data = df_all)
isSingular(model)
summary(model)


colnames(df_all)
lm_mat <- data.frame(matrix(nrow = 0, ncol = 6))
for (sample in unique(df_all$sample)) {
  print(sample)
  for (cluster in unique(df_all$clusters)) {
    sub_mat <-  df_all[df_all$sample == sample  & df_all$clusters == cluster, ]
    for (hallmark in paste0("H", 1:13)) {
      f <- formula(paste0(hallmark, " ~ score_scaled"))
      lm_out <- lm(f, sub_mat)
      lm_mat <- rbind(lm_mat, c(sample, hallmark, cluster, lm_out$coefficients[2],summary(lm_out)$adj.r.squared, mean(sub_mat[,hallmark])))
    }
  }
}


for (sample in unique(df_all$sample)) {
  print(sample)
  sub_mat <-  df_all[df_all$sample == sample, ]
  for (hallmark in paste0("H", 1:13)) {
    f <- formula(paste0(hallmark, " ~ score_scaled"))
    lm_out <- lm(f, sub_mat)
    lm_mat <- rbind(lm_mat, c(sample, hallmark, "all", lm_out$coefficients[2],summary(lm_out)$adj.r.squared))
  }
}




colnames(lm_mat) <- c("Sample", "Hallmark", "Cluster", "Slope", "Rsquared", "avg")

lm_barplot <- function(cluster, hallmark) {
  sub_lm <- lm_mat[lm_mat$Cluster==cluster & lm_mat$Hallmark == hallmark,]
  sub_lm$Rsquared2 <-  as.numeric(sub_lm$Rsquared)*as.numeric(sub_lm$Slope)/abs(as.numeric(sub_lm$Slope))

  tmp <- as.numeric(sub_lm$Rsquared2)
  names(tmp) <- sub_lm$Sample
  sub_lm$Sample <- factor(sub_lm$Sample, levels = names(sort(tmp)))
  return(ggplot(sub_lm, aes(x=Sample, y=as.numeric(Rsquared2), fill=as.numeric(avg))) + geom_bar(stat="identity") +
    ggtitle(paste0("Cluster: ",cluster, ", Hallmark: ", hallmark)) + theme_classic() + 
      scale_fill_gradient2(low = "blue", mid = "white",high = "red"))
}

h <- "H12"
lm_barplot(1, "H1")
lm_barplot(1, "H2")
lm_barplot(1, "H3")
lm_barplot(1, "H4")
lm_barplot(1, "H5")
lm_barplot(1, "H6")
lm_barplot(1, "H7")
lm_barplot(1, "H8")
lm_barplot(1, "H9")
lm_barplot(1, "H10")
lm_barplot(1, "H11")
lm_barplot(1, "H12")
lm_barplot(1, "H13")


lm_barplot(5, "H1")
lm_barplot(5, "H2")
lm_barplot(5, "H3")
lm_barplot(5, "H4")
lm_barplot(5, "H5")
lm_barplot(5, "H6")
lm_barplot(5, "H7")
lm_barplot(5, "H8")
lm_barplot(5, "H9")
lm_barplot(5, "H10")
lm_barplot(5, "H11")
lm_barplot(5, "H12")
lm_barplot(5, "H13")


lm_barplot(2, h)
lm_barplot(3, h)
lm_barplot(4, h)
lm_barplot(5, h)

lm_barplot("all", "H12")


ggplot(df_all, aes())



###### Model: Estiamte + score

lm_mat <- data.frame(matrix(nrow = 0, ncol = 16))
for (sample in unique(df_all$sample)) {
  print(sample)
  for (cluster in unique(df_all$clusters)) {
    sub_mat <-  df_all[df_all$sample == sample  & df_all$clusters == cluster, ]
    for (hallmark in paste0("H", 1:13)) {
      f <- formula(paste0(hallmark, " ~  estimate*score_scaled"))
      lm_out <- lm(f, sub_mat)
      lm_summary <- summary(lm_out)
      lm_mat <- rbind(lm_mat, c(sample, hallmark, cluster, lm_summary$adj.r.squared,
                                lm_summary$coefficients[2,1], lm_summary$coefficients[2,2], lm_summary$coefficients[2,3], lm_summary$coefficients[2,4],
                                lm_summary$coefficients[3,1], lm_summary$coefficients[3,2], lm_summary$coefficients[3,3], lm_summary$coefficients[3,4],
                                lm_summary$coefficients[4,1], lm_summary$coefficients[4,2], lm_summary$coefficients[4,3], lm_summary$coefficients[4,4]))
    }
  }
}





colnames(lm_mat) <- c("Sample", "Hallmark", "Cluster", "Rsquared", 
                      "EstimateSlope", "EstimateStdError", "EstimateTvalue", "EstimatePvalue",
                      "NeighborhoodSlope", "NeighborhoodStdError", "NeighborhoodTvalue", "NeighborhoodPvalue",
                      "InteractionSlope", "InteractionStdError", "InteractionTvalue", "InteractionPvalue")


hallmark <- "H12"
cluster <- "1"
sub_lm <- lm_mat[lm_mat$Cluster == cluster & lm_mat$Hallmark == hallmark,]
#sub_lm$Rsquared2 <-  as.numeric(sub_lm$Rsquared)*as.numeric(sub_lm$Slope)/abs(as.numeric(sub_lm$Slope))

tmp <- as.numeric(sub_lm$InteractionTvalue)
names(tmp) <- sub_lm$Sample
sub_lm$Sample <- factor(sub_lm$Sample, levels = names(sort(tmp)))
ggplot(sub_lm, aes(x=Sample, y=as.numeric(InteractionTvalue), fill=as.numeric(NeighborhoodSlope))) + geom_bar(stat="identity") +
  
  ggtitle(paste0("Cluster: ", cluster, ", Hallmark: ", hallmark)) + theme(axis.text.x = element_text(angle = 315, hjust = 0)) + 
    scale_fill_gradient2(low = "blue", mid = "white",high = "red")






lm_mat <- data.frame(matrix(nrow = 0, ncol = 16))
for (sample in unique(df_all$sample)) {
  print(sample)
  sub_mat <-  df_all[df_all$sample == sample  & df_all$clusters %in% c(1,5), ]
  sub_mat$clusters <- factor(sub_mat$clusters, levels = c(1,5))
  #sub_mat$clusters <- factor(sub_mat$clusters, levels = c(3,1,2,4,5))
  for (hallmark in paste0("H", 1:13)) {
    f <- formula(paste0(hallmark, " ~  score_scaled+clusters+clusters:score_scaled"))
    lm_out <- lm(f, sub_mat)
    lm_summary <- summary(lm_out)
    # lm_mat <- rbind(lm_mat, c(sample, hallmark, "all", lm_summary$adj.r.squared,
    #                           lm_summary$coefficients[2,1], lm_summary$coefficients[2,2], lm_summary$coefficients[2,3], lm_summary$coefficients[2,4],
    #                           lm_summary$coefficients[3,1], lm_summary$coefficients[3,2], lm_summary$coefficients[3,3], lm_summary$coefficients[3,4],
    #                           lm_summary$coefficients[4,1], lm_summary$coefficients[4,2], lm_summary$coefficients[4,3], lm_summary$coefficients[4,4],
    #                           lm_summary$coefficients[5,1], lm_summary$coefficients[5,2], lm_summary$coefficients[5,3], lm_summary$coefficients[5,4],
    #                           lm_summary$coefficients[6,1], lm_summary$coefficients[6,2], lm_summary$coefficients[6,3], lm_summary$coefficients[6,4],
    #                           lm_summary$coefficients[7,1], lm_summary$coefficients[7,2], lm_summary$coefficients[7,3], lm_summary$coefficients[7,4],
    #                           lm_summary$coefficients[8,1], lm_summary$coefficients[8,2], lm_summary$coefficients[8,3], lm_summary$coefficients[8,4],
    #                           lm_summary$coefficients[9,1], lm_summary$coefficients[9,2], lm_summary$coefficients[9,3], lm_summary$coefficients[9,4],
    #                           lm_summary$coefficients[10,1], lm_summary$coefficients[10,2], lm_summary$coefficients[10,3], lm_summary$coefficients[10,4]
    # )
    lm_mat <- rbind(lm_mat, c(sample, hallmark, "all", lm_summary$adj.r.squared,
                    lm_summary$coefficients[2,1], lm_summary$coefficients[2,2], lm_summary$coefficients[2,3], lm_summary$coefficients[2,4],
                    lm_summary$coefficients[3,1], lm_summary$coefficients[3,2], lm_summary$coefficients[3,3], lm_summary$coefficients[3,4],
                    lm_summary$coefficients[4,1], lm_summary$coefficients[4,2], lm_summary$coefficients[4,3], lm_summary$coefficients[4,4]
                    ))
  }
}

colnames(lm_mat) <- c("Sample", "Hallmark", "Cluster", "Rsquared", 
                      "NeighborhoodSlope", "NeighborhoodStdError", "NeighborhoodTvalue", "NeighborhoodPvalue",
                      "Cluster1Slope", "Cluster1StdError", "Cluster1Tvalue", "Cluster1Pvalue",
                      "Cluster2Slope", "Cluster2StdError", "Cluster2Tvalue", "Cluster2Pvalue",
                      "Cluster4Slope", "Cluster4StdError", "Cluster4Tvalue", "Cluster4Pvalue",
                      "Cluster5Slope", "Cluster5StdError", "Cluster5Tvalue", "Cluster5Pvalue",
                      "Interaction1Slope", "Interaction1StdError", "Interaction1Tvalue", "Interaction1Pvalue",
                      "Interaction2Slope", "Interaction2StdError", "Interaction2Tvalue", "Interaction2Pvalue",
                      "Interaction4Slope", "Interaction4StdError", "Interaction4Tvalue", "Interaction4Pvalue",
                      "Interaction5Slope", "Interaction5StdError", "Interaction5Tvalue", "Interaction5Pvalue")

colnames(lm_mat) <- c("Sample", "Hallmark", "Cluster", "Rsquared", 
                      "NeighborhoodSlope", "NeighborhoodStdError", "NeighborhoodTvalue", "NeighborhoodPvalue",
                      "Cluster5Slope", "Cluster5StdError", "Cluster5Tvalue", "Cluster5Pvalue",
                      "Interaction5Slope", "Interaction5StdError", "Interaction5Tvalue", "Interaction5Pvalue")

hallmark <- "H12"
sub_lm <- lm_mat[lm_mat$Hallmark == hallmark,]
#sub_lm$Rsquared2 <-  as.numeric(sub_lm$Rsquared)*as.numeric(sub_lm$Slope)/abs(as.numeric(sub_lm$Slope))

tmp <- as.numeric(sub_lm$InteractionTvalue)
names(tmp) <- sub_lm$Sample
sub_lm$Sample <- factor(sub_lm$Sample, levels = names(sort(tmp)))
ggplot(sub_lm, aes(x=Sample, y=as.numeric(InteractionTvalue), fill=as.numeric(NeighborhoodSlope))) + geom_bar(stat="identity") +
         
         ggtitle(paste0("Cluster: ","all", ", Hallmark: ", hallmark)) + theme(axis.text.x = element_text(angle = 315, hjust = 0)) + 
  scale_fill_gradient2(low = "blue", mid = "white",high = "red") +
  ggplot(sub_lm, aes(x=Sample, y=as.numeric(InteractionTvalue), fill=as.numeric(EstimateSlope))) + geom_bar(stat="identity") +
  
  ggtitle(paste0("Cluster: ","all", ", Hallmark: ", hallmark)) + theme(axis.text.x = element_text(angle = 315, hjust = 0)) + 
  scale_fill_gradient2(low = "blue", mid = "white",high = "red")



#DU8, UKF269T, P270. HCC5D
ggplot(df_all[df_all$sample=="Colorectal",], aes(x=score, y=H12)) + geom_smooth(method = "lm") + facet_wrap(~clusters) +
ggplot(df_all[df_all$sample=="DU8",], aes(x=estimate, y=H12)) + geom_smooth(method = "lm") + facet_wrap(~clusters)

ggplot(df_all[df_all$sample=="UKF269T",], aes(x=score, y=H12)) + geom_smooth(method = "lm") + facet_wrap(~clusters)
ggplot(df_all[df_all$sample=="UKF269T",], aes(x=estimate, y=H12)) + geom_smooth(method = "lm") + facet_wrap(~clusters)

ggplot(df_all[df_all$sample=="P270",], aes(x=score, y=H12)) + geom_smooth(method = "lm") + facet_wrap(~clusters) +
ggplot(df_all[df_all$sample=="P270",], aes(x=estimate, y=H12)) + geom_smooth(method = "lm") + facet_wrap(~clusters)

ggplot(df_all[df_all$sample=="ICC1L",], aes(x=score, y=H12)) + geom_smooth(method = "lm") + facet_wrap(~clusters)
ggplot(df_all[df_all$sample=="ICC1L",], aes(x=estimate, y=H12)) + geom_smooth(method = "lm") + facet_wrap(~clusters)


####################

#How many slopes are significant for each hallmark?
table(as.numeric(lm_mat$NeighborhoodPvalue) < 0.05, lm_mat$Hallmark)
lm_mat$signif <- as.numeric(as.numeric(lm_mat$NeighborhoodPvalue) < 0.05)

ggplot(lm_mat[as.numeric(lm_mat$NeighborhoodPvalue) < 0.05/13,], aes(x=Hallmark)) + geom_bar() + 
  geom_hline(yintercept = length(unique(lm_mat$Sample)), linetype="dashed") + theme_classic()

#How many slopes are significant for each sample?
ggplot(lm_mat[as.numeric(lm_mat$NeighborhoodPvalue) < 0.05/13/41,], aes(x=Sample)) + geom_bar() + geom_hline(yintercept = 13, linetype="dashed") + theme_classic() +
  theme(axis.text.x = element_text(angle = 315, hjust = 0)) 



#
hallmark <- "H12"
sub_lm <- lm_mat[lm_mat$Hallmark == hallmark,]
#sub_lm$Rsquared2 <-  as.numeric(sub_lm$Rsquared)*as.numeric(sub_lm$Slope)/abs(as.numeric(sub_lm$Slope))
tmp <- as.numeric(sub_lm$Interaction5Tvalue)
names(tmp) <- sub_lm$Sample
sub_lm$Sample <- factor(sub_lm$Sample, levels = names(sort(tmp)))
ggplot(sub_lm, aes(x=Sample, y=as.numeric(Interaction5Tvalue), fill=as.numeric(Cluster5Slope))) + geom_bar(stat="identity") +
  
  ggtitle(paste0("Cluster: ", cluster, ", Hallmark: ", hallmark)) + theme(axis.text.x = element_text(angle = 315, hjust = 0)) + 
  scale_fill_gradient2(low = "blue", mid = "white",high = "red")

