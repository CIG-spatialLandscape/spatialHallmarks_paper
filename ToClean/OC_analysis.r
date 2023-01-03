library(SeuratDisk)

OC <- readRDS("OC/OC_SC/OVSC_izar2020.st.comb.rds")


SaveH5Seurat(OC, "OC_log1p.h5Seurat", overwrite = T)
Convert("OC_log1p.h5Seurat", dest= "h5ad", overwrite = T)

OC <- RunPCA(OC, assay = "SCT", verbose = FALSE)
OC <- FindNeighbors(OC, reduction = "pca", dims = 1:50)
OC <- FindClusters(OC, verbose = FALSE)
OC <- RunUMAP(OC, reduction = "pca", dims = 1:50)



####################
library(Seurat)
library(patchwork)
library(ggplot2)

DimPlot(OC, reduction = "umap", group.by = "MajorLabels")
DimPlot(OC, reduction = "umap", group.by = "SubTypes")


####################
p = 'OC/spatial/Parent_Visium_Human_OvarianCancer_filtered_feature_bc_matrix.h5'

expr.data <- Seurat::Read10X_h5(filename = p)
OC.st <- Seurat::CreateSeuratObject(counts = expr.data, project = 'oc', assay = 'Spatial', na.rm=T)



img <- Seurat::Read10X_Image(image.dir = 'OC/spatial')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = OC.st)]
OC.st[['OC']] <- img

OC.st <- SCTransform(OC.st, assay = "Spatial", verbose = TRUE)

plot1 <- VlnPlot(OC.st, features = 'nCount_SCT', pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(OC.st, features = 'nCount_SCT') + theme(legend.position = "right")
wrap_plots(plot1, plot2)

plot1 <- VlnPlot(OC.st, features = 'nCount_Spatial', pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(OC.st, features = 'nCount_Spatial') + theme(legend.position = "right")
wrap_plots(plot1, plot2)

OC.st <- RunPCA(OC.st, assay = "SCT", verbose = FALSE)
OC.st <- FindNeighbors(OC.st, reduction = "pca", dims = 1:50)
OC.st <- FindClusters(OC.st, verbose = FALSE, resolution = 0.7)
OC.st <- RunUMAP(OC.st, reduction = "pca", dims = 1:50)


#OC.st <- FindSubCluster(OC.st, graph.name ="SCT_snn", resolution = 1.2)
colors = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
           "#44AA99", "#999933", "#882255", "#661100", "#6699CC")
p1 <- DimPlot(OC.st, reduction = 'umap', label = TRUE)
p2 <- SpatialDimPlot(OC.st, label = TRUE, label.size = 3, group.by = "seurat_clusters", cols = colors)
p1 + p2


#############
OC.st <- readRDS("OC/RDS/oc_st.rds")
names(OC.st[[]])


OC_sce <- readRDS("OC/RDS/OC_sce.rds")
OC_sce$col <- -OC_sce$col
palette <- RColorBrewer::brewer.pal(8, "Paired")
clusterPlot(OC_sce, palette = palette, size=0.05) + labs(title = "Spot-level clustering")
OC_sce$spatial.cluster
OC.st$bayes <- OC_sce$spatial.cluster[match(rownames(OC.st[[]]),colnames(OC_sce))]
SpatialDimPlot(OC.st, label = TRUE, label.size = 3, group.by = "bayes", cols = colors)



################
DefaultAssay(OC.st) <- "SCT"
OC.st <- NormalizeData(OC.st)
samples.combined <- SCTransform(samples.combined, assay = "RNA")
DefaultAssay(samples.combined) <- "SCT"
anchors <- FindTransferAnchors(reference = samples.combined, query = OC.st, normalization.method = "SCT", n.trees=1000, reduction = "cca")

predictions.assay <- TransferData(anchorset = anchors, refdata = samples.combined$MajorType, prediction.assay = TRUE,
                                  weight.reduction = "cca", dims = 1:30, n.trees = 1000)

OC.st[["prediction_majortypes"]] <- predictions.assay
DefaultAssay(OC.st) <- "prediction_majortypes"
p <- SpatialFeaturePlot(object = OC.st, features = rownames(OC.st), , ncol = 3)
pdf(file = "Desktop/IJC/datasets/IGTP/4A/figures/seurat_transfer_main.pdf")
plot(p)
dev.off()

samples.combined$SubType <- Idents(samples.combined)
samples.combined$MajorType <- as.character(Idents(samples.combined))
samples.combined$MajorType[samples.combined$SubType %in% c("CAFs_Immune1", "CAFs_Immune2")] <- "CAFs"
samples.combined$MajorType[samples.combined$SubType %in% c("Treg", "CD4+Tcells", "CD8+Tcells", "ProliferativeTcells")] <- ("Tcells"
DefaultAssay(OC.st) <- "SCT"
feat=markers[markers$cluster==4,]$gene[1:10]
SpatialFeaturePlot(OC.st, features = feat)


OC_sce <- readRDS("OC/RDS/OC_sce.rds")
OC_sce.enhanced <- spatialEnhance(OC_sce, q=8, d=15, platform="Visium", nrep=100000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE)

samples.combined <- FindSubCluster(samples.combined, graph.name = "integrated_snn", cluster = "MKI67+ ESCs", resolution = 0.1)
samples.combined$sub.cluster[samples.combined$sub.cluster %in% c("Epithelial_3", "MKI67+ ESCs_0")] <- "Proliferative"
samples.combined$sub.cluster[samples.combined$sub.cluster %in% c("Epithelial_0", "Epithelial_1", "Epithelial_2", "Epithelial_3")] <- "Malignant"
