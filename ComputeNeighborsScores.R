##################################################
## Project: Cancer Hallmarks
## Script purpose: Compute Neighborhood score for each sub-spot and create a data frame alongside hallmark activities
## Date: 22/12/2022
## Author: Sergi Cervilla * & Mustafa Sibai *
##################################################

## Auxiliary functions
# sort clusters from Cancer(1) to TME(5)
new_labels <- function(original, order) {
  labels <- original
  # Cancer (lowest)
  labels[original == order[1]] <- "1"
  labels[original == order[2]] <- "2"
  labels[original == order[3]] <- "3"
  labels[original == order[4]] <- "4"
  # TME (highest)
  labels[original == order[5]] <- "5"
  return(labels)
}
# select neighbours within a radius (thershold)
select_neighbours <- function(spot, threshold, distances) {
  return(grep(x = colnames(distances)[distances[spot, ] < threshold], pattern = spot, value = T, invert = T))
}

## CODE HERE
sample <- ""
# load hallmarks activities and exclude first sub-spot
hallmarks <- read.table("", sep = "\t", header = T)
hallmarks <- hallmarks[rownames(hallmarks) != "subspot_1.1", ]
# load coordinates activities and exclude first sub-spot
coord <- read.table("", sep = "\t", header = T)
coord <- coord[rownames(coord) != "subspot_1.1", ]
# load estimate values (no values for the first sub-spot)
estimate <- read.table("", sep = "\t", skip = 2, header = T)
estimate <- t(estimate[3, 3:ncol(estimate)])

source("../utils/CoordinatesEnhanced.R")
# compute real distances for sub-spots
coord <- subspot_dist(coord)

##### select neighbors
distances <- as.matrix(dist(coord[, c("realrow", "realcol")]))
gc()

# compute neighbors scores for each of the spots
score <- sapply(rownames(hallmarks), function(spot) {
  neigh <- select_neighbours(spot, 200, distances)
  score <- mean(estimate[neigh, ] / distances[spot, neigh])
})

# combine hallmark activities and neighbors scores
df <- data.frame(score = score)
df <- cbind(df, scale(hallmarks))
# split estimate values into 5 clusters and order them
labels <- kmeans(estimate, centers = 5)$cluster
labels_order <- order(sapply(1:5, function(x) {
  mean(estimate[labels == x])
}))
# add metadata to the data frame: estiamte cluster, estiamte values and sample
df$clusters <- new_labels(labels, labels_order)
df$estimate <- estimate
df$sample <- sample
#save data frame
write.table(df, "", sep = "\t")
