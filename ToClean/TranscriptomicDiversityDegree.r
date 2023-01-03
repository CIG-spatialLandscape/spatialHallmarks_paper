STobject <- readRDS("Desktop/IJC/datasets/IGTP/4A/RDS/OV4A.rds")
library(tictoc)
library(parallel)
hvg <- VariableFeatures(STobject, assay = "SCT")[1:2000]

correlations <- c()

counts <- as.matrix(STobject@assays$SCT@data[hvg,])


cl <- makeCluster (16)
clusterExport (cl, varlist = c("counts"))
tic("parallel")
a <- sapply(1:(l-1), function(i) {
    print(i)
    parSapply(cl, (i+1):l, function(j) {
      cor(counts[,i], counts[,j])
    })
})
toc()
stopCluster (cl)

correlations <- do.call(c, a)

MAD <- function(x) {
  return(median(abs(x - median(x))))
}

transcriptome_diversity_degree <- 1.4826*MAD(correlations)


l <- ncol(STobject)
cl <- makeCluster (16)
clusterExport (cl, varlist = c("counts", "l"))
tic("parallel")
a <- parSapply(cl, 1:(l-1), function(i) {
  print(i)
  sapply((i+1):l, function(j) {
    cor(counts[,i], counts[,j])
  })
})
toc()
stopCluster (cl)


files <- list.files(path="Desktop/IJC/datasets/Dutreneo/RDS", pattern="*.rds", full.names=TRUE, recursive=FALSE)
files <- files[c(2,4,5,6)]
tdd_c <- c()

tic("tdd")
for (file in files) {
  STobject <- readRDS(file)
  hvg <- VariableFeatures(STobject)[1:2000]
  hvg <- hvg[!is.na(hvg)]
  counts <- as.matrix(STobject@assays$SCT@data[hvg,])
  l <- ncol(STobject)
  cl <- makeCluster (16)
  clusterExport (cl, varlist = c("counts", "l"))
  correlations <- parSapply(cl, 1:(l-1), function(i) {
    print(i)
    sapply((i+1):l, function(j) {
      cor(counts[,i], counts[,j])
    })
  })
  
  correlations <- do.call(c, correlations)
  tdd_c <- c(tdd_c,1.4826*MAD(correlations))
}
toc()
