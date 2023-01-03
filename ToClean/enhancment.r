##
library(ggplot2)
library(BayesSpace)


set.seed(100)
### 4A
OV_4A <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/4A_sce.rds")

OV_4A.enhanced <- spatialEnhance(OV_4A, q=10, d=15, platform="Visium", init=OV_4A$spatial.cluster, 
                                 nrep=200000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE)

palette <- RColorBrewer::brewer.pal(10, name = "Paired")
p <- clusterPlot(OV_4A, palette = palette, size=0.05) + labs(title = "SubSpot-level clustering")
ggsave("Desktop/IJC/datasets/IGTP/figuresPaper/figures/BayesSpace/4A_enhanced.png", width = 7, height = 7)
gc()

saveRDS(OV_4A.enhanced, "Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/4A_enhanced.rds")
rm(OV_4A, OV_4A.enhanced)
gc()


### 5A
OV_5A <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/5A_sce.rds")
OV_5A <- spatialCluster(OV_5A, q=9, d=15, platform = "Visium", init.method = "kmeans", nrep = 50000, gamma = 3)

OV_5A$spatial.cluster <- as.factor(OV_5A$spatial.cluster)



OV_5A.enhanced <- spatialEnhance(OV_5A, q=9, d=15, platform="Visium", init=OV_5A$spatial.cluster, 
                                 nrep=200000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE)

palette <- RColorBrewer::brewer.pal(9, name = "Paired")
p <- clusterPlot(OV_5A, palette = palette, size=0.05) + labs(title = "SubSpot-level clustering")
ggsave("Desktop/IJC/datasets/IGTP/figuresPaper/figures/BayesSpace/5A_enhanced.png", width = 7, height = 7)
gc()

saveRDS(OV_5A.enhanced, "Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/5A_enhanced.rds")
rm(OV_5A, OV_5A.enhanced)
gc()



### 5B
OV_5B <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/5B_sce.rds")


OV_5B <- spatialCluster(OV_5B, q=8, d=15, platform = "Visium", init.method = "kmeans", nrep = 50000, gamma = 3)
OV_5B$spatial.cluster <- as.factor(OV_5B$spatial.cluster)
OV_5B.enhanced <- spatialEnhance(OV_5B, q=8, d=15, platform="Visium", init=OV_5B$spatial.cluster, 
                                 nrep=200000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE)

palette <- RColorBrewer::brewer.pal(10, name = "Paired")
p <- clusterPlot(OV_5B, palette = palette, size=0.05) + labs(title = "SubSpot-level clustering")
ggsave("Desktop/IJC/datasets/IGTP/figuresPaper/figures/BayesSpace/5B_enhanced.png", width = 7, height = 7)
gc()

saveRDS(OV_5B.enhanced, "Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/5B_enhanced.rds")
rm(OV_5B, OV_5B.enhanced)
gc()



### Colorectal
Colorectal <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/Colorectal_sce.rds")
Colorectal.enhanced <- spatialEnhance(Colorectal, q=7, d=15, platform="Visium", init=Colorectal$spatial.cluster, 
                                 nrep=200000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE)

palette <- RColorBrewer::brewer.pal(7, name = "Paired")
p <- clusterPlot(Colorectal.enhanced, palette = palette, size=0.05) + labs(title = "SubSpot-level clustering")
ggsave("Desktop/IJC/datasets/Public/Colorectal/figures/subcluster.png", width = 7, height = 7)
gc()

saveRDS(Colorectal.enhanced, "Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/Colorectal_enhanced.rds")
rm(Colorectal, Colorectal.enhanced)
gc()




### Breast
Breast <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/Breast_sce.rds")
Breast.enhanced <- spatialEnhance(Breast, q=8, d=15, platform="Visium", init=Breast$spatial.cluster, 
                                      nrep=200000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE)

palette <- RColorBrewer::brewer.pal(8, name = "Paired")
p <- clusterPlot(Breast.enhanced, palette = palette, size=0.05) + labs(title = "SubSpot-level clustering")
ggsave("Desktop/IJC/datasets/Public/Breast/figures/subcluster.png", width = 7, height = 7)
gc()

saveRDS(Breast.enhanced, "Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/Breast_enhanced.rds")
rm(Breast, Breast.enhanced)
gc()


### Glioblastoma
Glioblastoma <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/Glioblastoma_sce11.rds")
Glioblastoma.enhanced <- spatialEnhance(Glioblastoma, q=11, d=15, platform="Visium", init=Glioblastoma$spatial.cluster, 
                                      nrep=200000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE)

palette <- RColorBrewer::brewer.pal(11, name = "Paired")
p <- clusterPlot(Glioblastoma.enhanced, palette = palette, size=0.05) + labs(title = "SubSpot-level clustering") 
ggsave("Desktop/IJC/datasets/Public/Glioblastoma/figures/subcluster.png", width = 7, height = 7)
gc()

saveRDS(Glioblastoma.enhanced, "Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/Glioblastoma_enhanced.rds")
rm(Glioblastoma, Glioblastoma.enhanced)
gc()


### Ductal
Ductal <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/Ductal_sce.rds")
Ductal.enhanced <- spatialEnhance(Ductal, q=11, d=15, platform="Visium", init=Ductal$spatial.cluster, 
                                      nrep=200000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE)

palette <- RColorBrewer::brewer.pal(11, name = "Paired")
p <- clusterPlot(Ductal.enhanced, palette = palette, size=0.05) + labs(title = "SubSpot-level clustering") + scale_x_reverse()
ggsave("Desktop/IJC/datasets/Public/Ductal/figures/subcluster.png", width = 7, height = 7)
gc()

saveRDS(Ductal.enhanced, "Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/Ductal_enhanced.rds")
rm(Ductal, Ductal.enhanced)
gc()


### Ovarian
Ovarian <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/Ovarian_sce.rds")
Ovarian.enhanced <- spatialEnhance(Ovarian, q=8, d=15, platform="Visium", init=Ovarian$spatial.cluster, 
                                      nrep=200000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE)

palette <- RColorBrewer::brewer.pal(8, name = "Paired")
p <- clusterPlot(Ovarian.enhanced, palette = palette, size=0.05) + labs(title = "SubSpot-level clustering") + scale_x_reverse()
ggsave("Desktop/IJC/datasets/Public/Ovarian/figures/subcluster.png", width = 7, height = 7)
gc()

saveRDS(Ovarian.enhanced, "Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/Ovarian_enhanced.rds")
rm(Ovarian, Ovarian.enhanced)
gc()


### IC
IC <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/IC_sce.rds")
IC.enhanced <- spatialEnhance(IC, q=8, d=15, platform="Visium", init=IC$spatial.cluster, 
                                   nrep=200000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE)

palette <- RColorBrewer::brewer.pal(8, name = "Paired")
p <- clusterPlot(IC.enhanced, palette = palette, size=0.05) + labs(title = "SubSpot-level clustering")
ggsave("Desktop/IJC/datasets/Public/IC/figures/subcluster.png", width = 7, height = 7)
gc()

saveRDS(IC.enhanced, "Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/IC_enhanced.rds")
rm(IC, IC.enhanced)
gc()


