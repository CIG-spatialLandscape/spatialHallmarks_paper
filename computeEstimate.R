

library(Seurat)
library(estimate)

#read enhanced object with gene expression in the Seurat object format
STobject <- ""

#save expression matrix in a gct file
write.table(as.metrix(STobject@assays$SCT@data), "matrix.gct", quote = F, sep = "\t")

### Run bash script CreateHeaderGCT.sh

input = "matrix_input.gct"
output = "matrix_output.gct"
estimateScore(input, output)

#assign estiamte scores to each sub-spot in the seurat object
gct <- read.table(output, sep = "\t", header = F, skip = 3)
STobject@images[[1]] <- NULL
STobject <- STobject[,-1]


#cluster the estimate scores and sort them from cancer (1) to TME (5)
labels <- kmeans(t(gct[3,3:ncol(gct)]), centers = 5)$cluster

labels_order <- order(sapply(1:5, function(x){
  mean(t(gct[3,3:ncol(gct)])[labels==x])
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

STobject$estimate.cluster <- new_labels(labels,labels_order)


# Compute the average hallmark activity within each cluster (after running HallmarkScores.R)
avg_act <- matrix(nrow = 5, ncol = 13)
rownames(avg_act) <- paste0(name,"_", 1:5)
colnames(avg_act) <-  paste0("H",1:13)

for (cluster in 1:5) {
  for (hallmark in paste0("H",1:13)) {
    avg_act[cluster, hallmark] <- mean(STobject@meta.data[STobject$estimate.cluster == cluster, hallmark])
  }
}
#save average hallmark activities
write.table(avg_act,  "", quote = FALSE, sep = "\t", row.names = T)
#save Seurat object
saveRDS(STobject, "")

