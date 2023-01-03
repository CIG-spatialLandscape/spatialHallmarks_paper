###########################################

#         Olbercht 2021 dataset

##########################################
library(Seurat)
library(dplyr)
################################################################################
#                         Load data
################################################################################


#
counts <- readRDS("Desktop/IJC/datasets/olbrecht2021/RDS/2095-Olbrecht_counts_matrix.rds")
metadata <- read.csv("Desktop/IJC/datasets/olbrecht2021/RDS/2093-Olbrecht_metadata.csv")

samples_merged <- CreateSeuratObject(counts)
samples_merged$sample_name <- sapply(colnames(counts), function(x) substr(x, start = 18, stop = 25))
samples_merged@meta.data <-left_join(x = samples_merged@meta.data, y = metadata, by="sample_name")
rownames(samples_merged@meta.data) <- colnames(samples_merged)

################################################################################
#                         QC and filtering
################################################################################

mito.genes <- grep("^MT-", rownames(samples_merged@assays$RNA@data), value = TRUE)
print(paste(sort(mito.genes), collapse = ", "))

percent.mito <- Matrix::colSums(samples_merged@assays$RNA@counts[mito.genes, ]) /
  Matrix::colSums(samples_merged@assays$RNA@counts)
stopifnot(all.equal(names(percent.mito), colnames(samples_merged)))
samples_merged <- AddMetaData(object = samples_merged,
                 metadata = percent.mito,
                 col.name = "percent.mito")

selected_cells <- colnames(samples_merged)[samples_merged$percent.mito < 0.15]
samples_merged <- subset(samples_merged, cells = selected_cells)
selected_cells <- WhichCells(samples_merged, expression = nFeature_RNA < 6000)
samples_merged <- subset(samples_merged, cells = selected_cells)
selected_cells <- WhichCells(samples_merged, expression = nFeature_RNA > 200 )
samples_merged <- subset(samples_merged, cells = selected_cells)

selected_cells <- colnames(samples_merged)[samples_merged$sample_type =="Tumor"]
samples_merged <- subset(samples_merged, cells = selected_cells)

selected_cells <- colnames(samples_merged)[samples_merged$sample_site =="peritoneum"]
samples_merged <- subset(samples_merged, cells = selected_cells)

################################################################################
#                         Check batch effects
################################################################################
samples_merged <- NormalizeData(object = samples_merged,
                                     normalization.method = "LogNormalize",
                                     scale.factor = 10000)
samples_merged <- ScaleData(samples_merged)


samples_merged <- FindVariableFeatures(samples_merged)
samples_merged <- RunPCA(samples_merged, verbose = FALSE)


pca=samples_merged@reductions$pca
eival= (pca@stdev)^2
varexplaiend = eival/sum(eival)
plot(1:50, cumsum(varexplaiend), type="l")


samples_merged <- RunUMAP(samples_merged, reduction = "pca", dims = 1:30)

#clustering
samples_merged <- FindNeighbors(samples_merged, reduction = "pca", dims = 1:30)
samples_merged <- FindClusters(samples_merged, resolution = 0.35)
DimPlot(samples_merged, reduction = "umap")
DimPlot(samples_merged, reduction = "umap", group.by = "patient_id", repel=TRUE) +
DimPlot(samples_merged, reduction = "umap", group.by = "sample_type", repel=TRUE) +
DimPlot(samples_merged, reduction = "umap", group.by = "sample_site", repel=TRUE)
saveRDS(samples_merged, "Desktop/IJC/datasets/olbrecht2021/RDS/merged_tumorALL.rds")

################################################################################
#                         integration
################################################################################
samples_list <- SplitObject(samples_merged, split.by="patient_id")

samples_list <- lapply(X = samples_list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method="LogNormalize", scale.factor=10000)
  x <- ScaleData(x)
})

#Select variable features among all datasets
features <- SelectIntegrationFeatures(object.list = samples_list, nfeatures = 3000)


#select anchors from datasets
samples.anchors <- FindIntegrationAnchors(object.list = samples_list, normalization.method = "LogNormalize",
                                          anchor.features = features)
#rm(patient.list)
#saveRDS(patient.anchors, file="anchors.rds")

#integrate datasets
samples.combined <- IntegrateData(anchorset = samples.anchors, normalization.method = "LogNormalize")
DefaultAssay(samples.combined) <- "integrated"
#dimensional reduction
samples.combined <- ScaleData(samples.combined)
samples.combined <- RunPCA(samples.combined, verbose = FALSE, seed.use = NULL)

pca=samples.combined@reductions$pca
eival= (pca@stdev)^2
varexplaiend = eival/sum(eival)
plot(1:50, cumsum(varexplaiend), type="l")


samples.combined <- RunUMAP(samples.combined, reduction = "pca", dims = 1:50)

#clustering
samples.combined <- FindNeighbors(samples.combined, reduction = "pca", dims = 1:50)
samples.combined <- FindClusters(samples.combined, resolution =0.4)
seq(0,0.5, by=0.05)


DimPlot(samples.combined, reduction = "umap", repel=TRUE,label = T)

DimPlot(samples.combined, reduction = "umap", group.by = "patient_id", repel=TRUE) +
  DimPlot(samples.combined, reduction = "umap", group.by = "sample_type", repel=TRUE) +
  DimPlot(samples.combined, reduction = "umap", group.by = "sample_site", repel=TRUE)

saveRDS(samples.combined, "Desktop/IJC/datasets/olbrecht2021/RDS/integrated_peritoneum.rds")



#############
#all tumors
samples.combined <- RenameIdents(samples.combined, `0`="CAFs", `1`="Epithelial", `2`="Macrophages", `3`="Epithelial", 
                                 `4`="Fibroblast-P3", `5`="CAFs", `6`="Epithelial_Proliferative", `7`="T cells", 
                                 `8` = "Endothelial cells", `9` = "Fibroblast-Proliferative", `10` = "MyoFibroblast", `11`="Plasma B cells",
                                 `12` = "Epithelial", `13` = "CAFs", `14` = "CAFs", `15`="DC",
                                 `16` = "pDC", `17` = "Unassigned", `18` = "Endothelial cells")



################## Macrophages #####################################################
samples.combined <- FindSubCluster(samples.combined, cluster = "Macrophages", graph.name = "integrated_snn", resolution = 0.1)
DimPlot(samples.combined, reduction = "umap", label = T, group.by = "sub.cluster")
Idents(samples.combined) <- samples.combined$sub.cluster

m <- FindAllMarkers(samples.combined, only.pos = T)
top20 <- m %>% group_by(cluster) %>% top_n(20, avg_log2FC)
DoHeatmap(samples.combined, features = top20$gene)

sub <- subset(samples.combined, idents = c("Macrophages_0", "Macrophages_1", "Macrophages_2"))

samples.combined <- RenameIdents(samples.combined, `Macrophages_0`="TAMs",`Macrophages_1`="M1 Macrophages", `Macrophages_2`="TAMs")
################## TAMs #####################################################
samples.combined <- FindSubCluster(samples.combined, cluster = "T cells", graph.name = "integrated_snn", resolution = 0.2)
DimPlot(samples.combined, reduction = "umap", label = T)
Idents(samples.combined) <- samples.combined$sub.cluster

m <- FindAllMarkers(samples.combined, only.pos = T)
top20 <- m %>% group_by(cluster) %>% top_n(20, avg_log2FC)
DoHeatmap(samples.combined, features = top20$gene)

sub <- FindSubCluster(samples.combined, cluster = "T cells", graph.name = "integrated_snn", resolution = 0.2)
sub <- subset(samples.combined, idents = c("T cells"))
DimPlot(sub, reduction = "umap", label = T, group.by = "sub.cluster")



HGT_markers <- RunCellHGT(sub, pathways = markers, dims = 1:50, n.features = 200)
markers_gs_prediction <- rownames(HGT_markers)[apply(HGT_markers, 2, which.max)]
markers_gs_prediction_signif <- ifelse(apply(HGT_markers, 2, max)>2, yes = markers_gs_prediction, "unassigned")
sub$cellid_markers <- markers_gs_prediction_signif
DimPlot(sub, reduction = "umap", group.by = "cellid_markers", label = T)


ref.anchors <- FindTransferAnchors(reference = ref, query = sub,
                                   dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = ref.anchors, refdata = ref$cell_type,
                            dims = 1:30)
sub <- AddMetaData(sub, metadata = predictions)

DimPlot(sub, reduction = "umap", group.by = "predicted.id", label = T) + 
  FeaturePlot(sub, features = "prediction.score.max")

