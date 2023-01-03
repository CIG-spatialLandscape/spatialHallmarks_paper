##################################################
## Project: Cancer Hallmarks
## Script purpose: Compute TME Radar scores for Cancer spots
## Date: 22/12/2022
## Author: Sergi Cervilla & Mustafa Sibai
##################################################

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
x <- args[1] #sample


#load hallmark information for each spot
spot_info <- read.table("", sep = "\t")
#extract spot id to match with the coordinates
spot_info$spotid <- sapply(rownames(spot_info), function(spot){
  x <- strsplit(spot,  ".", fixed = T)
  paste0(x[[1]][1], ".", substr(x[[1]][2],1, 1))
})

#load spatial coordinates
coord <- read.table("", sep = "\t")
coord <- coord[spot_info$spotid,]
#add real coordinates to sub-spots
source("../utils/CoordinatesEnhanced.R")
coord <- subspot_coord(coord)
##
tme_spots <- spot_info$spotid[spot_info$clusters %in% c(4,5)]
cancer_spots <- spot_info$spotid[spot_info$clusters %in% c(1,2)]
#compute distances across TME and Cancer spots
coord <- coord[c(tme_spots, cancer_spots), c("realrow", "realcol")]
distances <- as.matrix(dist(coord))
distances <- distances[tme_spots, cancer_spots]
#create an empty data frame
tmp <- c()
#for each Cancer hallmark, apply the formula
for (hallmark in c("H1", "H3", "H5", "H6","H7", "H9", "H13")) {
  print(paste0(hallmark, "_", x))
  tmp <- rbind(tmp, sapply(cancer_spots, function(spot){
    sum(1/distances[, spot]*spot_info[spot_info$sample==x & spot_info$clusters %in% c(4,5), hallmark])
  }))
}
write.table(t(tmp),"", sep = "\t")