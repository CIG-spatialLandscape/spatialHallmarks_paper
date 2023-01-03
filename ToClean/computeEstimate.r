



ST_sce <- readRDS("Desktop/IJC/datasets/Public/Colorectal/RDS/Colorectal_sce.rds")  
ST_sce.enhanced <- spatialEnhance(ST_sce, q=7, d=15, platform="Visium", init=ST_sce$spatial.cluster, 
                                  nrep=20000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE)



saveRDS(ST_sce.enhanced, "Desktop/experiment/Colorectal_20000.rds")
rm(ST_sce, ST_sce.enhanced)
gc()


library(BayesSpace)
files <- list.files("Desktop/experiment", full.names = T)
for (file in files) {
  sce <- readRDS(file)
  p <- clusterPlot(sce)
  png(paste0(file, ".png"))
  plot(p)
  dev.off()
}



i <- 

SpatialFeaturePlot(STobject, features = paste0("H",i), pt.size.factor = 0.65) + scale_fill_gradientn("", colours = viridis::inferno(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 30)) + ggtitle(hallmark_names[i])


SpatialFeaturePlot(STobject, features = c("H1", "H_new1"), pt.size.factor = 0.65) & scale_fill_gradientn("", colours = viridis::inferno(100))
SpatialFeaturePlot(STobject, features = c("H2", "H_new2"), pt.size.factor = 0.65) & scale_fill_gradientn("", colours = viridis::inferno(100))
SpatialFeaturePlot(STobject, features = c("H3", "H_new3"), pt.size.factor = 0.65) & scale_fill_gradientn("", colours = viridis::inferno(100))
SpatialFeaturePlot(STobject, features = c("H4", "H_new4"), pt.size.factor = 0.65) & scale_fill_gradientn("", colours = viridis::inferno(100))
SpatialFeaturePlot(STobject, features = c("H5", "H_new5"), pt.size.factor = 0.65) & scale_fill_gradientn("", colours = viridis::inferno(100))
SpatialFeaturePlot(STobject, features = c("H6", "H_new6"), pt.size.factor = 0.65) & scale_fill_gradientn("", colours = viridis::inferno(100))
SpatialFeaturePlot(STobject, features = c("H7", "H_new7"), pt.size.factor = 0.65) & scale_fill_gradientn("", colours = viridis::inferno(100))
SpatialFeaturePlot(STobject, features = c("H8", "H_new8"), pt.size.factor = 0.65) & scale_fill_gradientn("", colours = viridis::inferno(100))
SpatialFeaturePlot(STobject, features = c("H9", "H_new9"), pt.size.factor = 0.65) & scale_fill_gradientn("", colours = viridis::inferno(100))
SpatialFeaturePlot(STobject, features = c("H10", "H_new11"), pt.size.factor = 0.65) & scale_fill_gradientn("", colours = viridis::inferno(100))
SpatialFeaturePlot(STobject, features = c("H11", "H_new11"), pt.size.factor = 0.65) & scale_fill_gradientn("", colours = viridis::inferno(100))
SpatialFeaturePlot(STobject, features = c("H12", "H_new12"), pt.size.factor = 0.65) & scale_fill_gradientn("", colours = viridis::inferno(100))
SpatialFeaturePlot(STobject, features = c("H13", "H_new13"), pt.size.factor = 0.65) & scale_fill_gradientn("", colours = viridis::inferno(100))





library(stringr)
files <- list.files("Desktop/IJC/datasets/IGTP/figuresPaper/distances_fixed", full.names = F)
files <- str_remove(files, "_distances.txt")

for(sample in files) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  write.table(assays(sce)$SCT, paste0("Desktop/estimate/input/input_", sample, ".gct"), sep = "\t", quote = F)
  rm(sce)
  gc()
}

library(estimate)
estimateScore("Desktop/estimate/input/tmp_input_Acinar.gct", output.ds = "Desktop/estimate/output/output_Acinar.gct")
estimateScore("Desktop/estimate/input/tmp_input_BreastA.gct", output.ds = "Desktop/estimate/output/output_BreastA.gct")
estimateScore("Desktop/estimate/input/tmp_input_Breast.gct", output.ds = "Desktop/estimate/output/output_Breast.gct")
estimateScore("Desktop/estimate/input/tmp_input_C20.gct", output.ds = "Desktop/estimate/output/output_C20.gct")
estimateScore("Desktop/estimate/input/tmp_input_C7.gct", output.ds = "Desktop/estimate/output/output_C7.gct")
estimateScore("Desktop/estimate/input/tmp_input_C34.gct", output.ds = "Desktop/estimate/output/output_C34.gct")
estimateScore("Desktop/estimate/input/tmp_input_C51.gct", output.ds = "Desktop/estimate/output/output_C51.gct")
estimateScore("Desktop/estimate/input/tmp_input_C21.gct", output.ds = "Desktop/estimate/output/output_C21.gct")
estimateScore("Desktop/estimate/input/tmp_input_cHC1T.gct", output.ds = "Desktop/estimate/output/output_cHC1T.gct")
estimateScore("Desktop/estimate/input/tmp_input_Colorectal.gct", output.ds = "Desktop/estimate/output/output_Colorectal.gct")
estimateScore("Desktop/estimate/input/tmp_input_CRC1.gct", output.ds = "Desktop/estimate/output/output_CRC1.gct")
estimateScore("Desktop/estimate/input/tmp_input_CRC2.gct", output.ds = "Desktop/estimate/output/output_CRC2.gct")
estimateScore("Desktop/estimate/input/tmp_input_CUP295.gct", output.ds = "Desktop/estimate/output/output_CUP295.gct")
estimateScore("Desktop/estimate/input/tmp_input_DU2.gct", output.ds = "Desktop/estimate/output/output_DU2.gct")
estimateScore("Desktop/estimate/input/tmp_input_DU3.gct", output.ds = "Desktop/estimate/output/output_DU3.gct")
estimateScore("Desktop/estimate/input/tmp_input_DU8.gct", output.ds = "Desktop/estimate/output/output_DU8.gct")
estimateScore("Desktop/estimate/input/tmp_input_DU12.gct", output.ds = "Desktop/estimate/output/output_DU12.gct")
estimateScore("Desktop/estimate/input/tmp_input_DU13.gct", output.ds = "Desktop/estimate/output/output_DU13.gct")
estimateScore("Desktop/estimate/input/tmp_input_Ductal.gct", output.ds = "Desktop/estimate/output/output_Ductal.gct")
estimateScore("Desktop/estimate/input/tmp_input_DuctalFFPE.gct", output.ds = "Desktop/estimate/output/output_DuctalFFPE.gct")
estimateScore("Desktop/estimate/input/tmp_input_Glioblastoma.gct", output.ds = "Desktop/estimate/output/output_Glioblastoma.gct")
estimateScore("Desktop/estimate/input/tmp_input_HCC1T.gct", output.ds = "Desktop/estimate/output/output_HCC1T.gct")
estimateScore("Desktop/estimate/input/tmp_input_HCC2T.gct", output.ds = "Desktop/estimate/output/output_HCC2T.gct")
estimateScore("Desktop/estimate/input/tmp_input_HCC5D.gct", output.ds = "Desktop/estimate/output/output_HCC5D.gct")
estimateScore("Desktop/estimate/input/tmp_input_IC.gct", output.ds = "Desktop/estimate/output/output_IC.gct")
estimateScore("Desktop/estimate/input/tmp_input_ICC1L.gct", output.ds = "Desktop/estimate/output/output_ICC1L.gct")
estimateScore("Desktop/estimate/input/tmp_input_Intestine.gct", output.ds = "Desktop/estimate/output/output_Intestine.gct")
estimateScore("Desktop/estimate/input/tmp_input_OV4A.gct", output.ds = "Desktop/estimate/output/output_OV4A.gct")
estimateScore("Desktop/estimate/input/tmp_input_Ovarian.gct", output.ds = "Desktop/estimate/output/output_Ovarian.gct")
estimateScore("Desktop/estimate/input/tmp_input_OVFFPE.gct", output.ds = "Desktop/estimate/output/output_OVFFPE.gct")
estimateScore("Desktop/estimate/input/tmp_input_OVD1.gct", output.ds = "Desktop/estimate/output/output_OVD1.gct")
estimateScore("Desktop/IJC/datasets/IGTP/figuresPaper/estimate/input/tmp_input_TNBCA.gct", output.ds = "Desktop/IJC/datasets/IGTP/figuresPaper/estimate/output/output_TNBCA.gct")

estimateScore("Desktop/estimate/input/tmp_input_UKF242T.gct", output.ds = "Desktop/estimate/output/output_UKF242T.gct")
estimateScore("Desktop/estimate/input/tmp_input_UKF260T.gct", output.ds = "Desktop/estimate/output/output_UKF260T.gct")
estimateScore("Desktop/estimate/input/tmp_input_UKF269T.gct", output.ds = "Desktop/estimate/output/output_UKF269T.gct")
estimateScore("Desktop/estimate/input/tmp_input_UKF275T.gct", output.ds = "Desktop/estimate/output/output_UKF275T.gct")

estimateScore("Desktop/estimate/input/tmp_input_P259_H2A2.gct", output.ds = "Desktop/estimate/output/output_P259_H2A2.gct")
estimateScore("Desktop/estimate/input/tmp_input_P264.gct", output.ds = "Desktop/estimate/output/output_P264.gct")
estimateScore("Desktop/estimate/input/tmp_input_P270.gct", output.ds = "Desktop/estimate/output/output_P270.gct")
estimateScore("Desktop/estimate/input/tmp_input_P288.gct", output.ds = "Desktop/estimate/output/output_P288.gct")
estimateScore("Desktop/estimate/input/tmp_input_P306.gct", output.ds = "Desktop/estimate/output/output_P306.gct")


estimateScore("Desktop/IJC/datasets/IGTP/figuresPaper/estimate/input/tmp_input_P3.gct", output.ds = "Desktop/IJC/datasets/IGTP/figuresPaper/estimate/output/output_P3.gct")
estimateScore("Desktop/IJC/datasets/IGTP/figuresPaper/estimate/input/tmp_input_P4.gct", output.ds = "Desktop/IJC/datasets/IGTP/figuresPaper/estimate/output/output_P4.gct")
estimateScore("Desktop/IJC/datasets/IGTP/figuresPaper/estimate/input/tmp_input_P6.gct", output.ds = "Desktop/IJC/datasets/IGTP/figuresPaper/estimate/output/output_P6.gct")
estimateScore("Desktop/IJC/datasets/IGTP/figuresPaper/estimate/input/tmp_input_PC1.gct", output.ds = "Desktop/IJC/datasets/IGTP/figuresPaper/estimate/output/output_PC1.gct")
estimateScore("Desktop/IJC/datasets/IGTP/figuresPaper/estimate/input/tmp_input_PC2.gct", output.ds = "Desktop/IJC/datasets/IGTP/figuresPaper/estimate/output/output_PC2.gct")
estimateScore("Desktop/IJC/datasets/IGTP/figuresPaper/estimate/input/tmp_input_Co2.gct", output.ds = "Desktop/IJC/datasets/IGTP/figuresPaper/estimate/output/output_Co2.gct")
estimateScore("Desktop/IJC/datasets/IGTP/figuresPaper/estimate/input/tmp_input_Co4.gct", output.ds = "Desktop/IJC/datasets/IGTP/figuresPaper/estimate/output/output_Co4.gct")
estimateScore("Desktop/IJC/datasets/IGTP/figuresPaper/estimate/input/tmp_input_M1.gct", output.ds = "Desktop/IJC/datasets/IGTP/figuresPaper/estimate/output/output_M1.gct")
estimateScore("Desktop/IJC/datasets/IGTP/figuresPaper/estimate/input/tmp_input_M2.gct", output.ds = "Desktop/IJC/datasets/IGTP/figuresPaper/estimate/output/output_M2.gct")
estimateScore("Desktop/IJC/datasets/IGTP/figuresPaper/estimate/input/tmp_input_M3.gct", output.ds = "Desktop/IJC/datasets/IGTP/figuresPaper/estimate/output/output_M3.gct")
estimateScore("Desktop/IJC/datasets/IGTP/figuresPaper/estimate/input/tmp_input_P1.gct", output.ds = "Desktop/IJC/datasets/IGTP/figuresPaper/estimate/output/output_P1.gct")
estimateScore("Desktop/IJC/datasets/IGTP/figuresPaper/estimate/input/tmp_input_P5.gct", output.ds = "Desktop/IJC/datasets/IGTP/figuresPaper/estimate/output/output_P5.gct")
estimateScore("Desktop/IJC/datasets/IGTP/figuresPaper/estimate/input/tmp_input_P8.gct", output.ds = "Desktop/IJC/datasets/IGTP/figuresPaper/estimate/output/output_P8.gct")
estimateScore("Desktop/IJC/datasets/IGTP/figuresPaper/estimate/input/tmp_input_Co1.gct", output.ds = "Desktop/IJC/datasets/IGTP/figuresPaper/estimate/output/output_Co1.gct")
estimateScore("Desktop/IJC/datasets/IGTP/figuresPaper/estimate/input/tmp_input_Co3.gct", output.ds = "Desktop/IJC/datasets/IGTP/figuresPaper/estimate/output/output_Co3.gct")
estimateScore("Desktop/IJC/datasets/IGTP/figuresPaper/estimate/input/tmp_input_M4.gct", output.ds = "Desktop/IJC/datasets/IGTP/figuresPaper/estimate/output/output_M4.gct")
estimateScore("Desktop/IJC/datasets/IGTP/figuresPaper/estimate/input/tmp_input_P7.gct", output.ds = "Desktop/IJC/datasets/IGTP/figuresPaper/estimate/output/output_P7.gct")


for(sample in files) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  write.table(as.matrix(GetAssayData(sce, assay = "SCT")), paste0("Desktop/estimate/input/input_", sample, ".gct"), sep = "\t", quote = F)
  rm(sce)
  gc()
}



for(sample in files) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  
  #spot label
  estimate <- read.table(paste0("Desktop/estimate/output/output_", sample, ".gct"), sep = "\t", skip = 2, header = T)
  
  sce <- sce[,-1]
  sce$estimate <- NA 
  sce$estimate <- t(estimate[3, 3:ncol(estimate)]) 
  sce$estimate_g <- kmeans(sce$estimate, centers = 2)$cluster 
  
  sce$estimate_g <- mclust::Mclust(sce$estimate, G = 3)$classification 
  #cluster label
  tmp <- read.table(paste0("Desktop/intersect/", sample, "/BS/estimate.gct"), sep = "\t", skip = 3, header = F)
  tmp <- data.frame(t(tmp[,-c(1:2)]), row.names = paste0(sample, "_", 1:(ncol(tmp)-2)))
  
  sce$cluster_score <- NA
  for (i in 1:10) {
    sce$cluster_score[sce$spatial.cluster == i] <- tmp[i, 3]
  }
  
  v <- .make_triangle_subspots(colData(sce), fill = "cluster_score")
  ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) + scale_fill_viridis_c(option="inferno") + theme_void()

  
  v <- .make_triangle_subspots(colData(sce), fill = "estimate_g")
  ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~factor(fill))) +  theme_void()

  v <- .make_triangle_subspots(colData(sce), fill = "estimate")
  ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~X3)) + scale_fill_viridis_c(option="inferno") + theme_void()
}

STobject$subspot.idx


set.seed(123)
sample <- "Colorectal"
estimate <- read.table(paste0("Desktop/estimate/output/output_", sample, ".gct"), sep = "\t", skip = 2, header = T)
estimate <- t(estimate[3, 3:ncol(estimate)])

labels <- kmeans(estimate, centers = 5)$cluster

labels_order <- order(sapply(1:5, function(x){
  mean(estimate[labels==x])
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

nl <- new_labels(labels,labels_order)

sce$estimate_g <- labels
v <- .make_triangle_subspots(colData(sce), fill = "estimate_g")
ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~factor(fill))) +  theme_void()

sce$estimate_g <- nl
v <- .make_triangle_subspots(colData(sce), fill = "estimate_g")
ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~factor(fill))) +  theme_void()
