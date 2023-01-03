library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

#### Create Seurat Object

setwd('Desktop/IJC/datasets/FFPE_PC_IF/')

p = 'FFPE_PC_IF/Visium_FFPE_Human_Prostate_IF_filtered_feature_bc_matrix.h5'

expr.data <- Seurat::Read10X_h5(filename = p)
prostate_if <- Seurat::CreateSeuratObject(counts = expr.data, project = 'prostate', assay = 'Spatial', na.rm=T)
prostate_if$slice <- 1
prostate_if$region <- 'prostate'

img.url = 'Visium_FFPE_Human_Prostate_IF_spatial.tar.gz'

img <- Seurat::Read10X_Image(image.dir = 'FFPE_PC_IF/spatial')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = prostate_if)]
prostate_if[['prostate']] <- img




#### QC - Pre-processing

plot1 <- VlnPlot(prostate_if, features = 'nCount_Spatial', pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(prostate_if, features = 'nCount_Spatial') + theme(legend.position = "right")
wrap_plots(plot1, plot2)

prostate_if <- SCTransform(prostate_if, assay = "Spatial", verbose = TRUE, return.only.var.genes = FALSE)

plot1_SCT <- VlnPlot(prostate_if, features = 'nCount_SCT', pt.size = 0.1) + NoLegend()
plot2_SCT <- SpatialFeaturePlot(prostate_if, features = 'nCount_SCT') + theme(legend.position = "right")
wrap_plots(plot1, plot2, plot1_SCT, plot2_SCT)


#### Dimension reduction

prostate_if <- RunPCA(prostate_if, assay = "SCT", verbose = FALSE)
prostate_if <- FindNeighbors(prostate_if, reduction = "pca", dims = 1:50)
prostate_if <- FindClusters(prostate_if, verbose = FALSE)
prostate_if <- RunUMAP(prostate_if, reduction = "pca", dims = 1:50)


prostate_if <- FindSubCluster(prostate_if, cluster = 7, graph.name ="SCT_snn", resolution = 1.2)
colors = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
           "#44AA99", "#999933", "#882255", "#661100", "#6699CC")

p1 <- DimPlot(prostate_if, reduction = 'umap', label = TRUE)
p2 <- SpatialDimPlot(prostate_if, label = TRUE, label.size = 3, group.by = "sub.cluster", )
p1 + p2
p2
prostate_if <- saveRDS(prostate_if, "prostate_seurat.rds")
prostate_if <- readRDS("prostate_seurat.rds")
prostate_sc <- readRDS("integrated_seurat.rds")

anchors <- FindTransferAnchors(reference = prostate_sc, query = prostate_if, normalization.method = "SCT", n.trees=1000, reduction = "cca", recompute.residuals = F)

rmpredictions.assay <- TransferData(anchorset = anchors, refdata = prostate_sc$sub.cluster, prediction.assay = TRUE,
                                  weight.reduction = "cca", dims = 1:30, n.trees = 1000)

prostate_if[["predictions_celltype"]] <- predictions.assay
SpatialFeaturePlot(object = prostate_if, features = c("Endothelial cells", "unassigned", "Fibroblasts", "Mast cells", "Basal cells", "Epithelial cells"), alpha = c(0.1, 1), ncol = 3)
SpatialFeaturePlot(object = prostate_if, features = c("Pericytes"), alpha = c(0.1, 1), ncol = 3)

DefaultAssay(prostate_if) <- "predictions_celltype"
p <- SpatialFeaturePlot(object = prostate_if, features = rownames(prostate_if), ncol = 3)
pdf(file = "seurat13_sub.pdf", width = 15, height = 40)
plot(p)
dev.off()
  
prostate_if <- FindSpatiallyVariableFeatures(prostate_if, assay = "SCT", features=VariableFeatures(prostate_if)[1:2000], selection.method="markvariogram")
features <- SpatiallyVariableFeatures(prostate_if, selection.method = "markvariogram")
