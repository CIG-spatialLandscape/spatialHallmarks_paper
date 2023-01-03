library(BayesSpace)
spatial_dir <- "Desktop/IJC/datasets/IGTP/4A/spatial/"


colData <- read.csv(file.path(spatial_dir, "tissue_positions_list.csv"), 
                    header = FALSE)
colnames(colData) <- c("spot", "in_tissue", "row", "col", 
                       "imagerow", "imagecol")
colData$spot <- paste0(colData$spot, "_1")
rownames(colData) <- colData$spot
colData <- colData[colData$in_tissue > 0, ]

### Convert to sce object
OV_4A <- SingleCellExperiment(assays =list(SCT = as.matrix(GetAssayData(OV4, slot = "data"))))
OV_4A <- SingleCellExperiment(assays =list(SCT = as.matrix(GetAssayData(OV4, slot = "data"))), colData = colData[colnames(OV_4A@assays@data$SCT),],)

metadata(OV_4A)$BayesSpace.data <- list()
metadata(OV_4A)$BayesSpace.data$platform <- "Visium"
metadata(OV_4A)$BayesSpace.data$is.enhanced <- FALSE

names(attributes(reducedDim(OV_4A)))
attributes(reducedDim(OV_4A))[[5]]
OV4A@reductions$NMF@cell.embeddings[colnames(OV_4A@assays@data$SCT),1:2]
reducedDim(OV_4A, "PCA")[,1:14] <-  OV4A@reductions$NMF@cell.embeddings[colnames(OV_4A@assays@data$SCT),1:14]


# Cluster with BayesSpace
set.seed(100)

OV_4A <- spatialPreprocess(OV_4A, platform = "Visium", n.PCs = 15,
                           n.HVGs = 2000, log.normalize = FALSE, assay.type = "SCT")

# Fin q for clusters
OV_4A <- qTune(OV_4A, qs = seq(2, 20), platform = "Visium",d =14)
p <- qPlot(OV_4A)

jpeg("Desktop/IJC/datasets/IGTP/4A/figures/BayesSpace/qPlot_filtered.jpeg", width = 8, height = 7, units = "in", res = 300)
print(p)
dev.off()

# 
q = 9
d = 15

Y <- reducedDim(OV_4A, "PCA")[, seq_len(d)]
set.seed(101)
init <- Mclust(Y[,-12], q, "EEE", verbose = FALSE)$classification

OV_4A$init <- init
OV_4A$init <- kmeans(Y, centers = 9)$cluster
# Run the clustering
set.seed(100)
OV_4A <- spatialCluster(OV_4A, q=q, d=d, platform = "Visium", nrep = 50000, gamma = 3)
  
OV_4A$spatial.cluster <- as.factor(OV_4A$spatial.cluster)

palette <- RColorBrewer::brewer.pal(q, name = "Paired")
clusterPlot(OV_4A, palette = palette, size=0.05) + labs(title = "Spot-level clustering")
clusterPlot(OV_4A, palette = palette, size=0.05, label = "init") + labs(title = "Spot-level clustering")
samples.combined


a <- OV_4A$spatial.cluster

OV_4A$clus <- OV4A$bayes.cluster2
OV4A$bayes.k <- OV_4A$spatial.cluster

OV4A  <- subset(OV4A, bayes.k %in% c(1,3,5,6,7,8))
