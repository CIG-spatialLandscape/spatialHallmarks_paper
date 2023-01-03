#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
sample <- args[1]

library(parallel)
library(dplyr)

###  Spatial continuity degree  ###
get_neighbors2 <- function(spot){
  row = coords[spot, "row"]
  col = coords[spot, "col"]
  return(c(rownames(coords[round(coords$row,7)==round(row-1/3,7) & round(coords$col,7)==round(col-1/3,7), ]),
           rownames(coords[round(coords$row,7)==round(row-1/3,7) & round(coords$col)==round(col+1/3,7), ]),
           rownames(coords[round(coords$row,7)==round(row,7) & round(coords$col)==round(col-2/3,7), ]),
           rownames(coords[round(coords$row,7)==round(row,7) & round(coords$col,7)==round(col+2/3,7), ]),
           rownames(coords[round(coords$row,7)==round(row+1/3,7) & round(coords$col,7)==round(col-1/3,7), ]),
           rownames(coords[round(coords$row,7)==round(row+1/3,7) & round(coords$col,7)==round(col+1/3,7), ])))
}

#extract data that will be used to compute
#spatial coordiantes
coords <-read.table(paste0("~/projects/Hallmarks/objects_mts/", sample, "_coords.txt"))
coords <- coords[rownames(coords) != "subspot_1.1",]
df <-  read.table(paste0("~/projects/Hallmarks/df_outs/", sample, ".txt"), sep = "\t")

df$spotid <- sapply(rownames(df), function(spot){
  x <- strsplit(spot,  ".", fixed = T)
  paste0(x[[1]][1], ".", substr(x[[1]][2],1, 1))
})
df$compartment <- "Cancer"
df$compartment[df$clusters %in% c(3)] <- "Buffer"
df$compartment[df$clusters %in% c(4,5)] <- "TME"
#cluster identities
clusters <- df$compartment
names(clusters) <- df$spotid

#sum of neighbours belonging to the same identity to the reference
n <- 0
#sum of all neighbours
d <- 0
#iterate by each spot of the tissue
for (spot in df$spotid[df$clusters != 3]) {
  #get neighbours id
  neigh <- get_neighbors2(spot)
  #iterate by each neighbour
  for (spot2 in neigh) {
    #add 1 if spot and its spot neighbour belong to the same cluster
    if (clusters[spot]==clusters[spot2]) n <- n+1
  }
  #add the number of neighbours by each spot
  d <- d + length(neigh)
}
#spatial continuity degree metric

write.table(n/d,paste0("~/projects/Hallmarks/SCD/", sample, ".txt"), sep = "\t")
