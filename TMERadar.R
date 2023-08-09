##################################################
## Project: Cancer Hallmarks
## Script purpose: Compute Radar scores of neoplastic Hallmarks in TME spots
## Author: Sergi Cervilla* & Mustafa Sibai*
##################################################

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

x <- args[1] #sample

#load hallmark information for each spot
spot_info <- read.table("", sep = "\t")

###### shift the data to positive numbers ####
subset_cols <- 1:13
df_subset <- spot_info[, subset_cols]

# Step 2: Find the minimum value in each column
min_vals <- sapply(df_subset, min)

# Step 3: Define a small constant value
constant <- 0.001

# Step 4: Subtract the minimum value and add the constant to each column
shifted_data <- data.frame(rownames(spot_info),
                           sapply(df_subset, function(x) x - min(x)) + constant,
                           stringsAsFactors = FALSE)

# Step 5: Assign column names to the shifted_data data.frame
names(shifted_data) <- c("rowname", colnames(df_subset))

# Step 6: Bind the remaining columns from the original data.frame to shifted_data
shifted_data <- cbind(shifted_data, spot_info[, -(subset_cols)])

rownames(shifted_data) <- shifted_data$rowname
shifted_data$rowname <- NULL
colnames(shifted_data)[14]<- "estimate.cluster"

spot_info <- shifted_data

######
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

tme_spots <- spot_info$spotid[spot_info$estimate.cluster %in% c(4,5)]
cancer_spots <- spot_info$spotid[spot_info$estimate.cluster %in% c(1,2)]
#compute distances across TME and Cancer spots
coord <- coord[c(tme_spots, cancer_spots), c("realrow", "realcol")]
distances <- as.matrix(dist(coord))
distances <- distances[tme_spots, cancer_spots]
#create an empty data frame
tmp <- c()
#for each Cancer hallmark, apply the formula
for (hallmark in paste0("H", c(2,4,8,9,10,11,12))) {
  print(paste0(hallmark, "_", x))
  tmp <- rbind(tmp, sapply(tme_spots, function(spot){
    sum(1/distances[spot,]*spot_info[spot_info$estimate.cluster %in% c(1,2), hallmark])
  }))
}

for  (hallmark in paste0("H", c(1,3,5,6,7,13))) {
  tmp <- rbind(tmp, sapply(tme_spots, function(spot){
    spot_info[spot_info$spotid == spot, hallmark]
  }))
}

tmp <- t(tmp)

colnames(tmp) <- c("H2_radar", "H4_radar", "H8_radar", "H9_radar", "H10_radar", "H11_radar", "H12_radar", "H1", "H3", "H5", "H6","H7", "H13")

write.table(tmp, "", sep = "\t")

