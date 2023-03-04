##################################################
## Project: Cancer Hallmarks
## Script purpose: Compute Spatial Continuity Degree at sub-spot resolution
## Date: 22/12/2022
## Author: Sergi Cervilla * & Mustafa Sibai *
##################################################

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
sample <- args[1]


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
#spatial coordinates to determine neighbors
coords <-read.table("")
coords <- coords[rownames(coords) != "subspot_1.1",]
#data frame containing sub-spot information
df <-  read.table("", sep = "\t")
#extract spot id to match with coordinates
df$spotid <- sapply(rownames(df), function(spot){
  x <- strsplit(spot,  ".", fixed = T)
  paste0(x[[1]][1], ".", substr(x[[1]][2],1, 1))
})

#Define Cancer, Buffer and TME compartments based on estimate clusters
df$compartment <- cut(df$clusters, breaks = c(0,2.5,3.5,5), labels = c("Cancer", "Buffer", "TME"))

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
#spatial continuity degree metric: n/d
write.table(n/d,"", sep = "\t")
