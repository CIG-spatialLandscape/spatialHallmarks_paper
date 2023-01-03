

library(SpatialPCA)

coord <- as.matrix(data.frame(x=ST@images$OC@coordinates$imagecol, y=-ST@images$OC@coordinates$imagerow, row.names = colnames(ST)))

ST_PCA = CreateSpatialPCAObject(counts=GetAssayData(ST, assay = "RNA", slot = "data"), location=coord, project = "SpatialPCA",
                            gene.type="spatial",sparkversion="spark", gene.number=3000,customGenelist=NULL,
                            min.loctions = 20, min.features=20, numCores_spark = 10)

ST_PCA = SpatialPCA_buildKernel(ST_PCA, kerneltype="gaussian", bandwidthtype="SJ")
ST_PCA = SpatialPCA_EstimateLoading(ST_PCA,fast=FALSE,SpatialPCnum=20)
ST_PCA = SpatialPCA_SpatialPCs(ST_PCA, fast=FALSE)

clusterlabel= walktrap_clustering(7, ST@SpatialPCs,round(sqrt(dim(ST@location)[1])))
clusterlabel_refine=refine_cluster_10x(clusterlabel,ST@location,shape="square")

# set color
cbp_spatialpca <- c(  "mediumaquamarine", "lightblue2","#F0E442",  "plum1","chocolate1","dodgerblue","palegreen4","red","#CC79A7","mediumpurple","seagreen1")
# visualize the cluster
plot_cluster(legend="right",location=ST@location,clusterlabel_refine,pointsize=5,text_size=20 ,title_in=paste0("SpatialPCA"),color_in=cbp_spatialpca)


library(slingshot)
# focus on tumor and its surrounding region
tumor_ind <- colnames(OC)[OC$bayes.cluster %in% c(1,2,5)]

sim<- SingleCellExperiment(assays = GetAssayData(OC, assay = "Spatial", slot = "data")[,tumor_ind])
reducedDims(sim) <- SimpleList(DRM = t(ST@SpatialPCs[,tumor_ind]))

colData(sim)$Walktrap <- factor(OC$bayes.cluster[OC$bayes.cluster %in% c(1,2,5)])    
colData(sim)$Walktrap <- factor(OC$bayes.cluster[OC$bayes.cluster %in% c(1,2,5)])  
sim  <-slingshot(sim, clusterLabels = 'Walktrap', reducedDim = 'DRM',start.clus="5", omega=T ) 
# in this data we set tumor region as start cluster
summary(sim@colData@listData)

pseudotime_traj1 = sim@colData@listData$slingPseudotime_1
gridnum = 10
cbp <- c(  "plum1", "chocolate1","dodgerblue","black")
p_traj1 = plot_trajectory(pseudotime_traj1, ST@location[tumor_ind,],OC$bayes.cluster[OC$bayes.cluster %in% c(1,2,5)],gridnum,cbp,pointsize=5 ,arrowlength=0.2,arrowsize=1,textsize=15 )

p_traj1$Arrowoverlay1
p_traj1$Pseudotime

STsimu_high_ST = SpatialPCA_highresolution(ST)
cluster_SpatialPCA_high = louvain_clustering(clusternum = 5, STsimu_high_ST, knearest = 3)
color_in=c(  "plum1", "palegreen4","mediumaquamarine",  "chocolate1","#F0E442","dodgerblue","lightblue2")
title_in="SpatialPCA High resolution"
plot_cluster(STsimu_high_ST@highPos, as.character(cluster_SpatialPCA_high), pointsize=2,text_size=20 ,title_in,color_in,legend="bottom")


oc_round = round(ST@location/100) # because we are using pixel coordinates provided by the original data paper, the scale of location coordinates is large, we divided by 100 for easily generating spots that cover the whole tissue slice.

x_in = c()
y_in = c()
for(y_coor in unique(loc_round[,2])){
  ind = which(loc_round[,2]==y_coor)
  x_ind = loc_round[ind,1]
  x_new = seq(from=min(x_ind),to=max(x_ind),by=1)
  y_new = rep(y_coor, length(x_new))
  x_in = c(x_in,x_new)
  y_in = c(y_in, y_new)
}
loc_round=data.frame(x_in,y_in)
x_in = c()
y_in = c()
for(x_coor in unique(loc_round[,1])){
  ind = which(loc_round[,1]==x_coor)
  y_ind = loc_round[ind,2]
  y_new = seq(from=min(y_ind),to=max(y_ind),by=1)
  x_new = rep(x_coor, length(y_new))
  y_in = c(y_in,y_new)
  x_in = c(x_in, x_new)
}
dat_sample = data.frame(x_in,y_in)*100

STsimu_high_ST_random = SpatialPCA_highresolution(ST,newloc = dat_sample)
cluster_SpatialPCA_high_custom = walktrap_clustering(8, latent_dat=STsimu_high_ST_random@highPCs,round(sqrt(dim(STsimu_high_ST_random@highPCs)[2])))
color_in=c(  "dodgerblue", "mediumaquamarine","#CC79A7",  "chocolate1","lightblue2","palegreen4","plum1","#F0E442","#CC79A7","mediumpurple","seagreen1")
title_in="SpatialPCA High resolution"
plot_cluster(STsimu_high_ST_random@highPos, as.character(cluster_SpatialPCA_high_custom), pointsize=1.5,text_size=20 ,title_in,color_in,legend="bottom")



#### Using PCA and 15 clusters
library(slingshot)
# focus on tumor and its surrounding region
tumor_ind <- colnames(OC)[OC$bayes.cluster %in% c(1,2,5)]

OV <- OV_4A[, tumor_ind]
#NMF <- OC@reductions$NMF@cell.embeddings[tumor_ind,]
#reducedDim(OV) <- NMF
OV  <-slingshot(OV, clusterLabels = 'spatial.cluster', reducedDim = 'PCA', start.clus = "14") 
# in this data we set tumor region as start cluster
summary(OV@colData@listData)

pseudotime_traj1 = OV@colData@listData$slingPseudotime_2
gridnum = 10
p_traj1 = plot_trajectory(pseudotime_traj1, ST@location[tumor_ind,],as.factor(OV$spatial.cluster),gridnum,pointsize=3.2 ,arrowlength=0.2,arrowsize=1,textsize=15, color_in = RColorBrewer::brewer.pal(11, name = "Set3"))

p_traj1$Arrowoverlay1
p_traj1$Pseudotime

OV@colData@listData$slingPseudotime_1


# focus on tumor and its surrounding region
tumor_ind <- colnames(OC)[OC$bayes.cluster %in% c(1,2,5)]

OV <- OV_4A[, tumor_ind]
OV$overlap <- OC[,tumor_ind]$experiment
#NMF <- OC@reductions$NMF@cell.embeddings[tumor_ind,]
#reducedDim(OV) <- NMF
OV  <-slingshot(OV, reducedDim = 'PCA', omega=T) 
# in this data we set tumor region as start cluster
summary(OV@colData@listData)

pseudotime_traj1 = OV@colData@listData$slingPseudotime_1
gridnum = 10
p_traj1 = plot_trajectory(pseudotime_traj1, ST@location[tumor_ind,],as.factor(OV$overlap),gridnum,pointsize=3.2 ,arrowlength=0.2,arrowsize=1,textsize=15, color_in = RColorBrewer::brewer.pal(11, name = "Set3"))

p_traj1$Arrowoverlay1
p_traj1$Pseudotime


OV$slingPseudotime_1



###########33

library(grDevices)
colors <- colorRampPalette(RColorBrewer::brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(OV$slingPseudotime_1, breaks=100)]

plot(reducedDims(OV)$PCA, col=OV$spatial.cluster,pch=16, asp = 1)
lines(SlingshotDataSet(OV), lwd=2, col='black')

ggplot(as.data.frame(reducedDims(OV)$PCA), aes(x=PC1, y=PC2, col=as.factor(OV$spatial.cluster))) + geom_point() + scale_color_manual(values = RColorBrewer::brewer.pal(11, name = "Set3")) +
lines(SlingshotDataSet(OV), lwd=2, col='black')


OC_subset <- OC[, tumor_ind]
phate <- phate(t(as.matrix(OC_subset@assays$SCT@data)), gamma = 1)


############3
reticulate::use_python("miniconda3/bin/python3", required = T)
reticulate::py_discover_config("phate")

ggplot(as.data.frame(phate$embedding), aes(x=PHATE1, y=PHATE2, col=as.factor(OV$spatial.cluster))) + geom_point() + scale_color_manual(values = RColorBrewer::brewer.pal(11, name = "Set3"))+ theme_classic() 
