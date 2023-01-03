###############
p <- qPlot(prostate_sce)
pdf(file="plot.pdf", width=15, height=40)
plot(p)
dev.off()
##################

library(BayesSpace)
library(ggplot2)
library(scuttle)

setwd("Desktop/IJC/datasets")

#read visium data
prostate_sce <- readVisium("FFPE_PC_IF")

#logNormalize
prostate_sce <- logNormCounts(prostate_sce)

#select top 2000 HVG
dec <- scran::modelGeneVar(prostate_sce)
top <- scran::getTopHVGs(dec, n=2000)

#preprocess
prostate_sce <- scater::runPCA(prostate_sce, subset_row=top)
prostate_sce <- spatialPreprocess(prostate_sce, platform = "Visium", skip.PCA = TRUE)

#qplot to estimate the number of clusters
prostate_sce <- qTune(prostate_sce, qs=seq(2,20))

q=10 #number of clusters (qplot)
d= 15 #number of PCA (15 as top performer)

#cluster at spot level
prostate_sce <- spatialCluster(prostate_sce, q=q, d=d, platform = "Visium", nrep = 50000, gamma = 3, save.chain = TRUE)
saveRDS(prostate_sce, "prostate_sce.rds")

#cluster at sub-spot level
prostate_sce.enhanced <- spatialEnhance(prostate_sce, q=q, d=d, platform="Visium", nrep=100000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE)
saveRDS(prostate_sce.enhanced, "prostate_sce.enhanced.rds")

palette <- RColorBrewer::brewer.pal(10, "Paired")
palette2 <- c(palette[8], palette[10],palette[3], palette[9], palette[2], palette[5], palette[7], palette[1], palette[6], palette[11], palette[4], palette[12])
#cluster plots
p1 <- clusterPlot(prostate_sce, palette = palette, size=0.05) + labs(title = "Spot-level clustering")
p2 <- clusterPlot(prostate_sce.enhanced, palette = palette, size=0.05) + labs(title = "SubSpot-level clustering") 
wrap_plots(p1, p2)

#enhance features
prostate_sce.enhanced <- enhanceFeatures(prostate_sce.enhanced, prostate_sce, feature_name=top, nrounds=0)

seurat_prostate <- as.Seurat(prostate_sce.enhanced, counts =NULL)
saveRDS(seurat_prostate, "seurat.enhanced.rds")

#########################################################################################################

library(dplyr)
library(stringr)
library(readr)

#load panglao markers db
panglao <- read_tsv("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz")

#filter and format the table
panglao <- panglao %>% filter(str_detect(species, "Hs"), `official gene symbol` %in% rownames(prostate_sce),
                                  `organ` %in% c("Connective tissue", "Epithelium", "Immune system", "Vasculature"))

panglao_v <- panglao %>% group_by(`cell type`) %>% summarise(geneset = list(`official gene symbol`))
all <- setNames(panglao_v$geneset, panglao_v$`cell type`)

#create a list storing markers for each cell type
markers <- list()
for (i in panglao_v$`cell type`) {
  markers[[i]] <- panglao_v$geneset[panglao_v$`cell type` == i][[1]]
}

#prostate_sce.enhanced <- enhanceFeatures(prostate_sce.enhanced, prostate_sce, model = "xgboost", feature_names = purrr::reduce(markers, c), nrounds=0)

sum_counts <- function(sce, features) {
  if (length(features) > 1) {
    colSums(logcounts(sce)[features, ])
  } else {
    logcounts(sce)[features, ]
  }
}
#sum the expression for each cell type
enhanced_expr <- purrr::map(markers, function(xs) sum_counts(prostate_sce.enhanced, xs))

plot_expression <- function(sce, expr, name) {
  featurePlot(sce, expr, color=NA) +
    viridis::scale_fill_viridis(option="A") +
    labs(title=name, fill="Log-normalized")
}

#plot for each cell type
enhanced_plots <- purrr::imap(enhanced_expr, function(x,y) plot_expression(prostate_sce.enhanced, x, y))
p <- patchwork::wrap_plots(enhanced_plots, ncol = 3)
pdf(file = "", width = 15, height = 40)
plot(p)
dev.off()

##################################################################################

seurat <- readRDS("prostate_seurat.rds")

sce<- as.SingleCellExperiment(seurat)

sce <- swapAltExp(sce, "SCT", saved = "original")
colData(sce) <- colData(sce1)
met

dec <- scran::modelGeneVar(sce)
top <- scran::getTopHVGs(dec, n=2000)
sce <- scater::runPCA(sce, subset_row=top)
sce <- spatialPreprocess(sce, platform = "Visium", log.normalize = FALSE, assay.type = "logcounts", n.HVGs = 2000, n.PCs = 15)
sce <- spatialCluster(sce, q=11, d=15, platform = "Visium", nrep = 50000, gamma = 3, save.chain = TRUE, init=sce@colData$spatial.cluster )

clusterPlot(sce, palette = palette, size=0.05) + labs(title = "Spot-level clustering")

sce@colData$spatial.cluster <- m[colnames(sce),]

sce <- SingleCellExperiment(assays=list(SCT=as.matrix(GetAssayData(prostate_if, slot="data"))))

m <- data.frame(row.names = colnames(sce1))
m <- data.frame(row.names = colnames(sce))
m$spatial.cluster <- sce@colData$seurat_clusters

m$row <- sce1@colData$row
m$col <- sce1@colData$col
m$imagerow <- sce1@colData$imagerow
m$imagecol <- sce1@colData$imagecol
m$sizeFactor <- sce1@colData$sizeFactor
m$spatial.cluster <- sce1@colData$spatial.cluster

sce@colData$row <- m[colnames(sce),1]
sce@colData$col <- m[colnames(sce),2]
sce@colData$imagerow <- m[colnames(sce),3]
sce@colData$imagecol <- m[colnames(sce),4]
sce@colData$sizeFactor <- m[colnames(sce),5]
sce@colData$spatial.clust <- m[colnames(sce),6]

palette2 <- c("")