#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(igraph)

source("/gpfs/scratch/bsc64/bsc64485/Hallmarks/Analysis/plotfunct.R")

annotation_compartments <- readRDS("/gpfs/scratch/bsc64/bsc64485/Hallmarks/Analysis/HPA/REAL_DISTANCE/annotation_compartments.rds")

df <- as.data.frame(matrix(nrow=0, ncol=5))

name <- args[1]

STobject <- readRDS(paste0("/gpfs/scratch/bsc64/bsc64485/Hallmarks/Analysis/k_factors/pipeline/objects/", name, "_enhanced_int.rds"))

STobject@meta.data[,paste0("H", 1:13)] <- scale(STobject@meta.data[,paste0("H", 1:13)])

if (name %in% c("UKF242T", "UKF260T", "UKF269T", "UKF275T"))  {
  coord <- STobject@images[[1]][[name]]@coordinates[, c("row", "col")] 
} else {coord <- STobject@images[[1]]@coordinates[, c("row", "col")]} 

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

sample = name
cancer <- c()
TME <- c()
Buffer <- c()
for (cluster in 1:length(unique(STobject$spatial.cluster))) {
 if (annotation_compartments[paste0(sample, "_", cluster)] == "Cancer") cancer <- c(cancer, colnames(STobject)[STobject$spatial.cluster == cluster])
 else if (annotation_compartments[paste0(sample, "_", cluster)] == "TME" ) TME <- c(TME, colnames(STobject)[STobject$spatial.cluster == cluster])
 else Buffer <- c(Buffer, colnames(STobject)[STobject$spatial.cluster == cluster])
}

#graph implementation
#cancer
adj_mat <- as.matrix(dist(coord[cancer, c("realrow", "realcol")], diag = T, upper = T))
adj_mat[adj_mat > 100] <- 0
g <- graph_from_adjacency_matrix(adj_mat, weighted=TRUE, mode = "undirected")
c <- components(g, mode = c("weak", "strong"))
component_id <- which(c$csize > 20)
cancer <- names(c$membership)[c$membership %in% component_id]
#tme
adj_mat <- as.matrix(dist(coord[TME, c("realrow", "realcol")], diag = T, upper = T))
adj_mat[adj_mat > 100] <- 0
g <- graph_from_adjacency_matrix(adj_mat, weighted=TRUE, mode = "undirected")
c <- components(g, mode = c("weak", "strong"))
component_id <- which(c$csize > 20)
TME <- names(c$membership)[c$membership %in% component_id]
#buffer
adj_mat <- as.matrix(dist(coord[Buffer, c("realrow", "realcol")], diag = T, upper = T))
adj_mat[adj_mat > 100] <- 0
g <- graph_from_adjacency_matrix(adj_mat, weighted=TRUE, mode = "undirected")
c <- components(g, mode = c("weak", "strong"))
component_id <- which(c$csize > 20)
Buffer <- names(c$membership)[c$membership %in% component_id]

#################################

#extract coordiantes
if (length(TME) != 0 & length(cancer) != 0) {
  distances <- as.matrix(dist(coord, diag = T, upper = T))
  distances <- distances[cancer, TME]
  
  ## Cancer
  min_dist <- apply(distances,1,min)
  for (hallmark in paste0("H", 1:13)) {
    tmp <- data.frame(d = min_dist, h=STobject@meta.data[cancer,hallmark])
    tmp$hallmark <- hallmark
    tmp$sample <- name
    tmp$h_type <- "Cancer_TME"
    df <- rbind(df, tmp)
  }
  ## TME
  min_dist <- apply(distances,2,min)
  for (hallmark in paste0("H", 1:13)) {
    tmp <- data.frame(d = min_dist, h=STobject@meta.data[TME,hallmark])
    tmp$hallmark <- hallmark
    tmp$sample <- name
    tmp$h_type <- "TME_Cancer"
    df <- rbind(df, tmp)
  }
}

#extract coordiantes
if (length(Buffer) != 0 & length(cancer) != 0) {
  distances <- as.matrix(dist(coord, diag = T, upper = T))
  distances <- distances[cancer, Buffer]
  
  ## Cancer
  min_dist <- apply(distances,1,min)
  for (hallmark in paste0("H", 1:13)) {
    tmp <- data.frame(d = min_dist, h=STobject@meta.data[cancer,hallmark])
    tmp$hallmark <- hallmark
    tmp$sample <- name
    tmp$h_type <- "Cancer_Buffer"
    df <- rbind(df, tmp)
  }
  ## Buffer
  min_dist <- apply(distances,2,min)
  for (hallmark in paste0("H", 1:13)) {
    tmp <- data.frame(d = min_dist, h=STobject@meta.data[Buffer,hallmark])
    tmp$hallmark <- hallmark
    tmp$sample <- name
    tmp$h_type <- "Buffer_Cancer"
    df <- rbind(df, tmp)
  }
}

#extract coordiantes
if (length(TME) != 0 & length(Buffer) != 0) {
  distances <- as.matrix(dist(coord, diag = T, upper = T))
  distances <- distances[Buffer, TME]
  
  ## Buffer
  min_dist <- apply(distances,1,min)
  for (hallmark in paste0("H", 1:13)) {
    tmp <- data.frame(d = min_dist, h=STobject@meta.data[Buffer,hallmark])
    tmp$hallmark <- hallmark
    tmp$sample <- name
    tmp$h_type <- "Buffer_TME"
    df <- rbind(df, tmp)
  }
  ## TME
  min_dist <- apply(distances,2,min)
  for (hallmark in paste0("H", 1:13)) {
    tmp <- data.frame(d = min_dist, h=STobject@meta.data[TME,hallmark])
    tmp$hallmark <- hallmark
    tmp$sample <- name
    tmp$h_type <- "TME_Buffer"
    df <- rbind(df, tmp)
  }
}

#rm(STobject)


write.table(df, paste0("/gpfs/scratch/bsc64/bsc64485/Hallmarks/Analysis/HPA/REAL_DISTANCE/", name, "_distances.txt"), sep = "\t")


