library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

#### Create Seurat Object

setwd('Desktop/IJC/datasets/FFPE_PC/')

p = 'FFPE_PC/Visium_FFPE_Human_Prostate_Cancer_filtered_feature_bc_matrix.h5'

expr.data <- Seurat::Read10X_h5(filename = p)
pc <- Seurat::CreateSeuratObject(counts = expr.data, project = 'prostate', assay = 'Spatial', na.rm=T)
pc$slice <- 1
pc$region <- 'prostate'

img <- Seurat::Read10X_Image(image.dir = 'FFPE_PC/spatial')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = pc)]
pc[['prostate']] <- img




#### QC - Pre-processing

plot1 <- VlnPlot(pc, features = 'nCount_Spatial', pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(pc, features = 'nCount_Spatial') + theme(legend.position = "right")
wrap_plots(plot1, plot2)

pc <- SCTransform(pc, assay = "Spatial", verbose = TRUE, return.only.var.genes = FALSE)

plot1_SCT <- VlnPlot(pc, features = 'nCount_SCT', pt.size = 0.1) + NoLegend()
plot2_SCT <- SpatialFeaturePlot(pc, features = 'nCount_SCT') + theme(legend.position = "right")
wrap_plots(plot1, plot2, plot1_SCT, plot2_SCT)


#### Dimension reduction

pc <- RunPCA(pc, assay = "SCT", verbose = FALSE)
pc <- FindNeighbors(pc, reduction = "pca", dims = 1:50)
pc <- FindClusters(pc, verbose = FALSE)
pc <- RunUMAP(pc, reduction = "pca", dims = 1:50)

colors = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
           "#44AA99", "#999933", "#882255", "#661100", "#6699CC")

p1 <- DimPlot(pc, reduction = 'umap', label = TRUE)
p2 <- SpatialDimPlot(pc, label = TRUE, label.size = 3, group.by = "seurat_clusters")
p1 + p2

saveRDS(pc, "pc_seurat.rds")
pc <- readRDS("prostate_seurat.rds")
prostate_sc <- readRDS("integrated_seurat.rds")

anchors <- FindTransferAnchors(reference = prostate_sc, query = pc, normalization.method = "SCT", n.trees=1000, reduction = "cca")

predictions.assay <- TransferData(anchorset = anchors, refdata = prostate_sc$sub.cluster, prediction.assay = TRUE,
                                  weight.reduction = "cca", dims = 1:30, n.trees = 1000)

pc[["predictions_celltype"]] <- predictions.assay
SpatialFeaturePlot(object = pc, features = c("Endothelial cells", "unassigned", "Fibroblasts", "Mast cells", "Basal cells", "Epithelial cells"), alpha = c(0.1, 1), ncol = 3)
SpatialFeaturePlot(object = pc, features = c("Pericytes"), alpha = c(0.1, 1), ncol = 3)

DefaultAssay(pc) <- "predictions_celltype"
p <- SpatialFeaturePlot(object = pc, features = rownames(pc), ncol = 3)
pdf(file = "seurat13_pc_sub.pdf", width = 15, height = 40)
plot(p)
dev.off()
