
library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)

source("Desktop/IJC/datasets/IGTP/figuresPaper/scripts/utilities/plotfunct.r")
source("Desktop/IJC/datasets/IGTP/figuresPaper/scripts/utilities/AnnotateCorrelations.r")

df <- as.data.frame(matrix(nrow=0, ncol=5))

files <- list.files(path="Desktop/enhanced", pattern="*.rds", full.names=TRUE, recursive=FALSE)

for (file in files) {
  STobject <- readRDS(file)
  print(file)
  name <- strsplit(strsplit(file, "/")[[1]][3], "_")[[1]][1] 
  if (length(strsplit(strsplit(file, "/")[[1]][3], "_")[[1]])==3) name <- paste0(name, "_FFPE")
  
  #real distances
  coord <- STobject@images[[1]]@coordinates[, c("row", "col")]
  
  
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
    }
    else if (n == 2) {
      realrow <- realrow + factor
      realcol <- realcol - factor
    }
    else if (n == 3) {
      realrow <- realrow - factor
      realcol <- realcol + factor
    }
    else if (n == 4) {
      realrow <- realrow - factor
      realcol <- realcol - factor
    }
    else if (n == 5) {
      realrow <- realrow 
      realcol <- realcol + 2 * factor
    }
    else if (n == 6) {
      realrow <- realrow 
      realcol <- realcol - 2 * factor
    }
    coord[subspot, c("realrow", "realcol")] <- c(realrow, realcol)
    
    
  }
  STobject@images[[1]]@coordinates[, c("realrow", "realcol")] <- coord[, c("realrow", "realcol")]
  
  
  
  factor <- STobject@meta.data[,paste0("factor_", 1:5)]
  #extract spots with highest factor activity for each factor
  factor_spots <- lapply(1:5, function(x) {
    colnames(STobject)[factor[,x] > mean(factor[,x]) + 1*sd(factor[,x])]
  })
  cancer_factors <- which(annotation_proximity[[name]] == "Cancer" )
  TME_factors <-which(annotation_proximity[[name]] == "TME")
  
  
  cancer <- sapply(cancer_factors, function(x){
    factor_spots[[x]]
  })
  if (class(cancer) == "list") cancer <- do.call(c, cancer)
  cancer <- cancer[!duplicated(cancer)]
  
  TME <- sapply(TME_factors, function(x){
    factor_spots[[x]]
  })
  if (class(TME) == "list") TME <- do.call(c, TME)
  TME <- TME[!duplicated(TME)]
  
  #remove overlapped spots
  duplicated <- c(TME, cancer)[duplicated(c(TME, cancer))]
  TME <- TME[!TME %in% duplicated]
  cancer <- cancer[!cancer %in% duplicated]
  
  #extract coordiantes
  coord <- STobject@images[[1]]@coordinates[c(cancer,TME),c("realrow", "realcol")]
  distances <- as.matrix(dist(coord, diag = T, upper = T))
  distances <- distances[cancer, TME]
  
  ## Cancer
  min_dist <- apply(distances,1,min)
  for (hallmark in paste0("H", 1:13)) {
    tmp <- data.frame(d = min_dist, h=scale(STobject@meta.data[cancer,hallmark]))
    tmp$hallmark <- hallmark
    tmp$sample <- name
    tmp$h_type <- "Cancer"
    df <- rbind(df, tmp)
  }
  ## TME
  min_dist <- apply(distances,2,min)
  for (hallmark in paste0("H", 1:13)) {
    tmp <- data.frame(d = min_dist, h=scale(STobject@meta.data[TME,hallmark]))
    tmp$hallmark <- hallmark
    tmp$sample <- name
    tmp$h_type <- "TME"
    df <- rbind(df, tmp)
  }
  #rm(STobject)
} 

write.table(df, "Desktop/...", sep = "\t")
  
  

#fitting
fitting <- as.data.frame(matrix(nrow=0, ncol=8))

for (tumor in unique(df$sample)) {
  for (hallmark in unique(df$hallmark)) {
    for (compartment in unique(df$h_type)) {
      df1 <- df[df$sample==tumor & df$hallmark==hallmark & df$h_type==compartment,]
      
      model <- loess(h ~ d, data=df1)
      xrange <- range(df1$d)
      xseq <- seq(from=xrange[1], to=xrange[2], length=500)
      pred <- predict(model, newdata = data.frame(d = xseq), se=T)
      y = pred$fit
      ci <- pred$se.fit * qt(0.95 / 2 + .5, pred$df)
      ymin = y - ci
      ymax = y + ci
      loess.DF <- data.frame(x = xseq, y, ymin, ymax, se = pred$se.fit)
      loess.DF$sample <- tumor
      loess.DF$hallmark <- hallmark
      loess.DF$compartment <- compartment
      fitting <- rbind(fitting, loess.DF)
      
    }
  }
}
  
write.table(fitting, "Desktop/...", sep = "\t")