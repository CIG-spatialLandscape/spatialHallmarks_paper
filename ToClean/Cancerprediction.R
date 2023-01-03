#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
x <- args[1]

df_all <- read.table(paste0("projects/Hallamrks/df_outs/", x, ".txt"), sep = "\t")



df_all$spotid <- sapply(rownames(df_all), function(spot){
  x <- strsplit(spot,  ".", fixed = T)
  paste0(x[[1]][1], ".", substr(x[[1]][2],1, 1))
})


df_cancer <- df_all[df_all$clusters %in% c(4,5),]

print(x)
coord <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/objects_mts/", x, "_coords.txt"), sep = "\t")
coord <- coord[df_all$spotid,]
#put real distances
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
##
tme_spots <- df_all$spotid[df_all$clusters %in% c(4,5)]
cancer_spots <- df_all$spotid[df_all$clusters %in% c(1,2)]
coord <- coord[c(tme_spots, cancer_spots), c("realrow", "realcol")]
distances <- as.matrix(dist(coord))
distances <- distances[tme_spots, cancer_spots]
gc()
tmp <- c()
for (hallmark in c("H1", "H3", "H5", "H6","H7", "H9", "H13")) {
  print(paste0(hallmark, "_", x))
  tmp <- rbind(tmp, sapply(cancer_spots, function(spot){
    sum(1/distances[,spot]*df_all[tme_spots, hallmark])
  }))
}

write.table(t(tmp),paste0("projects/Hallmarks/Cancer_predictors/", x, ".txt"), sep = "\t")




##


files <- list.files("Desktop/CancerPredictors/", full.names = T)
tmp <- lapply(files, function(x) {
  read.table(x, sep = "\t")

})
tmp <- do.call(rbind, tmp)
colnames(tmp) <- c("H1_TME", "H3_TME", "H5_TME", "H6_TME","H7_TME", "H9_TME", "H13_TME")
df_cancer <- df_all[df_all$clusters %in% c(1,2),]
df_cancer <- cbind(df_cancer, tmp)



df_cancer[is.na(df_cancer)] <- 0


