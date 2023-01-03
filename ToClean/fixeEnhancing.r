STobject <- readRDS("Desktop/enhanced/new/CUP_295_enhanced.rds")


STobject <- readRDS("Desktop/IJC/datasets/Dutreneo/RDS/Seurat_enhancement/DU2_enhanced.rds")
STobject <- FindVariableFeatures(STobject)
STobject <- RunPCA(STobject, assay = "SCT", verbose = FALSE)
STobject <- FindNeighbors(STobject, reduction = "pca", dims = 1:30)
STobject <- FindClusters(STobject, verbose = FALSE, resolution = 0.08)
SpatialDimPlot(STobject, pt.size.factor = 0.5)
SpatialDimPlot(STobject, pt.size.factor = 0.5, group.by = c("spatial.cluster", "seurat_clusters"))
VlnPlot(STobject, features = paste0("H",1:13))
VlnPlot(STobject, features = paste0("H",1:13), group.by = "spatial.cluster")
STobject <- RunUMAP(brain, reduction = "pca", dims = 1:30)



sce <- readRDS("Desktop/IJC/datasets/Dutreneo/RDS/sce/DU12_sce.rds")
clusterPlot(sce)

SpatialDimPlot(sce, group.by = "spatial.cluster", pt.size.factor = 0.5)




samples <- grep(list.files("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final", full.names = F,), pattern='*enhanced.rds', invert=F, value=TRUE)
samples <- grep(list.files("Desktop/IJC/datasets/Dutreneo/RDS/enhancement", full.names = F,), pattern='*.rds', invert=F, value=TRUE)
samples <- grep(list.files("Desktop/RDS_new", full.names = F,), pattern='*.rds', invert=F, value=TRUE)


palette <- DiscretePalette(15)
for (sample in samples) {
  sce <- readRDS(paste0("Desktop/RDS_new/", sample))
  png(paste0("Desktop/BS/", sample, ".png"))
  p <- clusterPlot(sce, palette = palette)
  plot(p)
  dev.off()
}










ST_sce <- readRDS(paste0("Desktop/IJC/datasets/Dutreneo/RDS/sce/DU12_sce.rds"))

ST_sce.enhanced <- spatialEnhance(ST_sce, q=length(unique(ST_sce$spatial.cluster)), d=15, platform="Visium", init=ST_sce$spatial.cluster, 
                                  nrep=1000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE, burn.in = 10)
clusterPlot(ST_sce.enhanced)


sce2 = ST_sce
sce2$imagecol.new = sce2$imagerow
sce2$imagerow = sce2$imagecol
sce2$imagecol = sce2$imagecol.new

ST_sce.enhanced2 <- spatialEnhance(sce2, q=length(unique(sce2$spatial.cluster)), d=15, platform="Visium", init=sce2$spatial.cluster, 
                                  nrep=20000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE, burn.in = 1000)
clusterPlot(ST_sce.enhanced2)

plot(ST_sce$row, ST_sce$imagerow)
plot(ST_sce$col, ST_sce$imagecol)


cor(ST_sce$col, ST_sce$imagerow)
cor(ST_sce$col, ST_sce$imagecol)








#DU12
ST_sce <- readRDS(paste0("Desktop/IJC/datasets/Dutreneo/RDS/sce/DU12_sce.rds"))
sce2 = ST_sce
sce2$imagecol.new = sce2$imagerow
sce2$imagerow = sce2$imagecol
sce2$imagecol = sce2$imagecol.new

ST_sce.enhanced <- spatialEnhance(sce2, q=length(unique(sce2$spatial.cluster)), d=15, platform="Visium", init=sce2$spatial.cluster, 
                                   nrep=20000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE, burn.in = 1000)
saveRDS(ST_sce.enhanced2, "Desktop/RDS_new/DU12_enhanced.rds")
rm(ST_sce, sce2, ST_sce.enhanced)




#DU2
ST_sce <- readRDS(paste0("Desktop/IJC/datasets/Dutreneo/RDS/sce/DU2_sce.rds"))
sce2 = ST_sce
sce2$imagecol.new = sce2$imagerow
sce2$imagerow = sce2$imagecol
sce2$imagecol = sce2$imagecol.new

ST_sce.enhanced <- spatialEnhance(sce2, q=length(unique(sce2$spatial.cluster)), d=15, platform="Visium", init=sce2$spatial.cluster, 
                                  nrep=200000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE, burn.in = 1000)
saveRDS(ST_sce.enhanced, "Desktop/RDS_new/DU2_enhanced.rds")
rm(ST_sce, sce2, ST_sce.enhanced)

#DU3
ST_sce <- readRDS(paste0("Desktop/IJC/datasets/Dutreneo/RDS/sce/DU3_sce.rds"))
sce2 = ST_sce
sce2$imagecol.new = sce2$imagerow
sce2$imagerow = sce2$imagecol
sce2$imagecol = sce2$imagecol.new

ST_sce.enhanced <- spatialEnhance(sce2, q=length(unique(sce2$spatial.cluster)), d=15, platform="Visium", init=sce2$spatial.cluster, 
                                  nrep=200000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE, burn.in = 1000)
saveRDS(ST_sce.enhanced, "Desktop/RDS_new/DU3_enhanced.rds")
rm(ST_sce, sce2, ST_sce.enhanced)

#DU8
ST_sce <- readRDS(paste0("Desktop/IJC/datasets/Dutreneo/RDS/sce/DU8_sce.rds"))
sce2 = ST_sce
sce2$imagecol.new = sce2$imagerow
sce2$imagerow = sce2$imagecol
sce2$imagecol = sce2$imagecol.new

ST_sce.enhanced <- spatialEnhance(sce2, q=length(unique(sce2$spatial.cluster)), d=15, platform="Visium", init=sce2$spatial.cluster, 
                                  nrep=200000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE, burn.in = 1000)
saveRDS(ST_sce.enhanced, "Desktop/RDS_new/DU8_enhanced.rds")
rm(ST_sce, sce2, ST_sce.enhanced)




#DU13
ST_sce <- readRDS(paste0("Desktop/IJC/datasets/Dutreneo/RDS/sce/DU13_sce.rds"))
sce2 = ST_sce
sce2$imagecol.new = sce2$imagerow
sce2$imagerow = sce2$imagecol
sce2$imagecol = sce2$imagecol.new

ST_sce.enhanced <- spatialEnhance(sce2, q=length(unique(sce2$spatial.cluster)), d=15, platform="Visium", init=sce2$spatial.cluster, 
                                  nrep=200000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE, burn.in = 1000)
saveRDS(ST_sce.enhanced, "Desktop/RDS_new/DU13_enhanced.rds")
rm(ST_sce, sce2, ST_sce.enhanced)

#ffpe_c_7_
ST_sce <- readRDS(paste0("Desktop/IJC/datasets/Public/ffpe_c_7/RDS/ffpe_c_7_sce.rds"))
sce2 = ST_sce
sce2$imagecol.new = sce2$imagerow
sce2$imagerow = sce2$imagecol
sce2$imagecol = sce2$imagecol.new

ST_sce.enhanced <- spatialEnhance(sce2, q=length(unique(sce2$spatial.cluster)), d=15, platform="Visium", init=sce2$spatial.cluster, 
                                  nrep=200000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE, burn.in = 1000)
saveRDS(ST_sce.enhanced, "Desktop/RDS_new/ffpe_c_7_enhanced.rds")
rm(ST_sce, sce2, ST_sce.enhanced)


#ffpe_c_20_
ST_sce <- readRDS(paste0("Desktop/IJC/datasets/Public/ffpe_c_20/RDS/ffpe_c_20_sce.rds"))
sce2 = ST_sce
sce2$imagecol.new = sce2$imagerow
sce2$imagerow = sce2$imagecol
sce2$imagecol = sce2$imagecol.new

ST_sce.enhanced <- spatialEnhance(sce2, q=length(unique(sce2$spatial.cluster)), d=15, platform="Visium", init=sce2$spatial.cluster, 
                                  nrep=200000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE, burn.in = 1000)
saveRDS(ST_sce.enhanced, "Desktop/RDS_new/ffpe_c_20_enhanced.rds")
rm(ST_sce, sce2, ST_sce.enhanced)

#ffpe_c_21_
ST_sce <- readRDS(paste0("Desktop/IJC/datasets/Public/ffpe_c_21/RDS/ffpe_c_21_sce.rds"))
sce2 = ST_sce
sce2$imagecol.new = sce2$imagerow
sce2$imagerow = sce2$imagecol
sce2$imagecol = sce2$imagecol.new

ST_sce.enhanced <- spatialEnhance(sce2, q=length(unique(sce2$spatial.cluster)), d=15, platform="Visium", init=sce2$spatial.cluster, 
                                  nrep=200000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE, burn.in = 1000)
saveRDS(ST_sce.enhanced, "Desktop/RDS_new/ffpe_c_21_enhanced.rds")
rm(ST_sce, sce2, ST_sce.enhanced)

#ffpe_c_34_
ST_sce <- readRDS(paste0("Desktop/IJC/datasets/Public/ffpe_c_34/RDS/ffpe_c_34_sce.rds"))
sce2 = ST_sce
sce2$imagecol.new = sce2$imagerow
sce2$imagerow = sce2$imagecol
sce2$imagecol = sce2$imagecol.new

ST_sce.enhanced <- spatialEnhance(sce2, q=length(unique(sce2$spatial.cluster)), d=15, platform="Visium", init=sce2$spatial.cluster, 
                                  nrep=200000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE, burn.in = 1000)
saveRDS(ST_sce.enhanced, "Desktop/RDS_new/ffpe_c_34_enhanced.rds")
rm(ST_sce, sce2, ST_sce.enhanced)

#ffpe_c_51
ST_sce <- readRDS(paste0("Desktop/IJC/datasets/Public/ffpe_c_51/RDS/ffpe_c_51_sce.rds"))
sce2 = ST_sce
sce2$imagecol.new = sce2$imagerow
sce2$imagerow = sce2$imagecol
sce2$imagecol = sce2$imagecol.new

ST_sce.enhanced <- spatialEnhance(sce2, q=length(unique(sce2$spatial.cluster)), d=15, platform="Visium", init=sce2$spatial.cluster, 
                                  nrep=200000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE, burn.in = 1000)
saveRDS(ST_sce.enhanced, "Desktop/RDS_new/ffpe_c_51_enhanced.rds")
rm(ST_sce, sce2, ST_sce.enhanced)


#CUP295
ST_sce <- readRDS(paste0("Desktop/IJC/datasets/Public/CUP295/RDS/CUP295_sce.rds"))
sce2 = ST_sce
sce2$imagecol.new = sce2$imagerow
sce2$imagerow = sce2$imagecol
sce2$imagecol = sce2$imagecol.new

ST_sce.enhanced <- spatialEnhance(sce2, q=length(unique(sce2$spatial.cluster)), d=15, platform="Visium", init=sce2$spatial.cluster, 
                                  nrep=200000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE, burn.in = 1000)
saveRDS(ST_sce.enhanced, "Desktop/RDS_new/CUP295_enhanced.rds")
rm(ST_sce, sce2, ST_sce.enhanced)