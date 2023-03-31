##################################################
## Project: Cancer Hallmarks
## Script purpose: Compute Cancer Radar scores for TME spots
## Author: Sergi Cervilla* & Mustafa Sibai*
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
for (hallmark in c("H2", "H4", "H8", "H10", "H11", "H12")) {
  print(paste0(hallmark, "_", x))
  tmp <- rbind(tmp, sapply(tme_spots, function(spot){
    sum(1/distances[spot,]*spot_info[spot_info$sample==x & spot_info$clusters %in% c(1,2), hallmark])
  }))
}
write.table(t(tmp), "", sep = "\t")

