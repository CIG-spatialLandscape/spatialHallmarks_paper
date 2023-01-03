###########################################

#         Izer 2020 dataset

##########################################

setwd("Desktop/IJC/datasets/")

################################################################################
#                         Load data
################################################################################


#

counts <- read.csv("OC/OC_SC/GSE146026_Izar_HGSOC_ascites_10x_log.tsv.gz", sep = "\t", row.names = 1)
OC.sc <- Seurat::CreateSeuratObject(counts = counts[8:ncol(counts),])
OC.sc$patient <- paste0("p",as.character(counts[2,]))
OC.sc$time <- as.integer(counts[3,])
OC.sc$sample_ID <- as.character(counts[4,])
OC.sc$clst <- as.integer(counts[5,])
rownames(OC.sc)


################################################################################
#                         QC and filtering
################################################################################
samples_list <- SplitObject(OC.sc, split.by = "patient")

for (i in 1:length(samples_list)) {
  samples_list[[i]] <- PercentageFeatureSet(samples_list[[i]], "^MT-", col.name = "percent_mito")
  samples_list[[i]] <- PercentageFeatureSet(samples_list[[i]], "^RP[SL]", col.name = "percent_ribo")

}

for (i in 1:length(samples_list)) {
  samples_list[[i]] <- PercentageFeatureSet(samples_list[[i]], "^MT.", col.name = "percent_mito")
  samples_list[[i]] <- PercentageFeatureSet(samples_list[[i]], "^RP[SL]", col.name = "percent_ribo")
  
  selected_cells <- WhichCells(samples_list[[i]], expression = nFeature_RNA > 600)
  samples_list[[i]] <- subset(samples_list[[i]], cells = selected_cells)
  total_counts_per_cell <- colSums(samples_list[[i]]@assays$RNA@counts)
  mt_genes <- rownames(samples_list[[i]])[grep("^MT.", rownames(samples_list[[i]]))]
  samples_list[[i]] <- samples_list[[i]][!grepl("^MT.", rownames(samples_list[[i]])),]
  samples_list[[i]] <- samples_list[[i]][!grepl("^RP[SL]", rownames(samples_list[[i]])),]
  
  samples_list[[i]][["RNA"]]@data <- expm1(samples_list[[i]][["RNA"]]@data)
  }

################################################################################
#                         Check batch effects
################################################################################


samples_merged <- merge(samples_list[[1]], y = samples_list[2:6], add.cell.ids = seq(1:length(samples_list)), project = "patient")

samples_merged <- SCTransform(samples_merged)

samples_merged <- RunPCA(samples_merged, verbose = FALSE)
samples_merged <- RunUMAP(samples_merged, reduction = "pca", dims = 1:30)

#clustering
samples_merged <- FindNeighbors(samples_merged, reduction = "pca", dims = 1:30)
samples_merged <- FindClusters(samples_merged, resolution = 0.3)

samples_merged$orig.ident[samples_merged$orig.ident == "X10x"] <- samples_merged$patient[samples_merged$orig.ident == "X10x"]
DimPlot(samples_merged, reduction = "umap", group.by = "orig.ident", repel=TRUE) + 
DimPlot(samples_merged, reduction = "umap", repel=TRUE)
saveRDS(samples_merged, "OC/RDS/merged_all.rds")

################################################################################
#                         integration
################################################################################

samples_list <- lapply(X = samples_list, FUN = SCTransform) 
#Select variable features among all datasets
features <- SelectIntegrationFeatures(object.list = samples_list, nfeatures = 3000)

samples_list <- PrepSCTIntegration(object.list = samples_list, anchor.features = features)

#select anchors from datasets
samples.anchors <- FindIntegrationAnchors(object.list = samples_list, normalization.method = "LogNormalize",
                                          anchor.features = features, reduction = "cca")
#rm(patient.list)
#saveRDS(patient.anchors, file="anchors.rds")

#integrate datasets
samples.combined <- IntegrateData(anchorset = samples.anchors, normalization.method = "LogNormalize")
DefaultAssay(samples.combined) <- "integrated"
samples.combined <- ScaleData(samples.combined)
#dimensional reduction
samples.combined <- RunPCA(samples.combined, verbose = FALSE, seed.use = NULL)
samples.combined <- RunUMAP(samples.combined, reduction = "pca", dims = 1:35)

#clustering
samples.combined <- FindNeighbors(samples.combined, reduction = "pca", dims = 1:35)
samples.combined <- FindClusters(samples.combined, resolution = 0.6)

p1 <- DimPlot(samples.combined, reduction = "umap", group.by = "patient") + labs(title = "Patient distribution")
p2 <- DimPlot(samples.combined, reduction = "umap", repel = TRUE, label = T) + labs(title = "Clusters (resolution=0.6)", caption = "data:Izar2020")
p1 + p2 

clustree(samples.combined, prefix = "integrated_snn_res.", node_colour="sc3_stability")

#saveRDS(samples.combined, "olalekan2021/RDS/integrated.rds")


################################################################################
#                         Manual annotation
################################################################################

markers <- FindAllMarkers(samples.combined, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
conserved_markers1 <- FindConservedMarkers(samples.combined, ident.1 = 1, grouping.var = "orig.ident")
samples.combined$orig.ident[samples.combined$orig.ident == "X10x"] <- samples.combined$patient[samples.combined$orig.ident == "X10x"]

feat = c("ACTA2", "CD44", "CDH2", "FN1")
feat = c("DCN", "COL1A1", "COL3A1", "CD4")
feat = c("CD4", "CD8", "IFNG", "FOXP3")
feat=rownames(markers)[1:10]
feat=rownames(conserved_markers1)[1:10]
feat=head(markers[markers$cluster==2,]$gene,20)
feat=markers[markers$cluster==5,]$gene[20:40]
VlnPlot(samples.combined, features = feat, ncol=3)
DotPlot(samples.combined, features = feat, group.by = "orig.ident")
DotPlot(samples.combined, features = feat)
FeaturePlot(samples.combined, features=feat)

feat = c("PDPN", "DCN", "THY1")
feat = c("CD1C", "CD1E", "CCR7", "CD83")
feat = c("CD19", "CD79A", "CD79B")
feat = c("GATA1")
### Annotation
# 0: Cancer related - need subclusters
# 1: Fibroblast - need subcluster
# 2: T cells
# 3: Macrophages - need sublsuter (bottom macrophages, top MHCII )
# 4: B cells
# 5: ? B cells 2 ?
# 6: Endothelial cells
# 7: ?
# 8: ?

major <- c("Cancer-Related", "Fibroblast", "T cells", "Macrophages", "B cells", "B cells (Plasma?)", "Endothelial cells", "U1", "U2")

samples.combined <- FindSubCluster(samples.combined, cluster = "Cancer-related", graph.name = "integrated_snn", resolution = 0.2)
DimPlot(samples.combined, reduction = "umap", repel = TRUE, label = T, group.by = "sub.cluster")
Idents(samples.combined) <- samples.combined$sub.cluster
samples.combined <- RenameIdents(samples.combined, `0` = "Cancer-related", `1` = "Fibroblast", `2` = "T cells", `3` = "Macrophages", `4` = "B cells", `5` = "B cells (plasma?)", `6` = "Endothelial cells",
                                 `7` = "U1", `8` = "U2")


markers_1 <- FindMarkers(samples.combined, ident.1 = "Cancer-related_0", ident.2 = c("Cancer-related_1", "Cancer-related_2"), only.pos = T)
markers_2 <- FindMarkers(samples.combined, ident.1 = "Cancer-related_1", ident.2 = c("Cancer-related_0", "Cancer-related_2"), only.pos = T)
markers_3 <- FindMarkers(samples.combined, ident.1 = "Cancer-related_2", ident.2 = c("Cancer-related_0", "Cancer-related_1"), only.pos = T)
feat = rownames(markers_1)[1:20]
feat = rownames(markers_2)[1:20]
feat = rownames(markers_3)[1:20]
markers1 <- FindAllMarkers(samples.combined, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
feat=head(markers1[markers1$cluster=="Cancer-related_0",]$gene,20)
feat=head(markers1[markers1$cluster=="Cancer-related_1",]$gene,20)
feat=head(markers1[markers1$cluster=="Cancer-related_2",]$gene,20)
feat=head(markers1[markers1$cluster=="Fibroblast_0",]$gene,20)
FeaturePlot(samples.combined, features=c("MUC16", "PAX8", "MKI67", "VIM", "ZEB1", "FN1"))
FeaturePlot(samples.combined, features=c("IFNGR1", "CD36", "DDX5", "MNDA", "C3AR1"))
FeaturePlot(samples.combined, features=c("EPCAM", "CD24", "CLDN3", "CLDN4", "CLDN7", "KRT8", "KRT19"))
FeaturePlot(samples.combined, features=c("CD74", "HLA.DRA", "HLA.DRB1", "HLA.DMA", "HLA.DQB1", "HLA.DRB5", "HLA.DMB", "HLA.DQA1", "HLA.DOA", "HLA.DQA2")) #MH2
DotPlot(samples.combined, features=c("CD74", "HLA-DRA", "HLA-DRB1", "HLA-DMA", "HLA-DQB1", "HLA-DRB5", "HLA-DMB", "HLA-DQA1", "HLA-DOA", "HLA-DQA2")) #MH2
DotPlot(samples.combined, features=c("CD14", "CD68", "PTPRC", "CD163", "CD53", "AIF1", "TYROBP"))

FeaturePlot(samples.combined, features=c("CD70", "CD90", "CD105", "CD34")) #MSCs


#####
samples.combined <- FindSubCluster(samples.combined, cluster = "Macrophages", graph.name = "integrated_snn", resolution = 0.1)
DimPlot(samples.combined, reduction = "umap", group.by = "sub.cluster", label = T)
samples.combined$sub.cluster[samples.combined$sub.cluster=="Macrophages_2"] <- "Macrophages_1"
feat=head(markers1[markers1$cluster=="Fibroblast_0",]$gene,20)


epithelial <- c("KRT17", "KRT6A", "KLK10", "KLK7", "KLK8", "KRT4", "EPCAM", "MMP7", "SOX17", "CXCL17", "CLDN7")
macro <- c("CD14", "S100A9", "S100A8", "C1QB", "C1QA", "AIF1", "FCER1G", "CCL3", "FCGR1A", "CSF1R")

CAF <- c("FGF5", "CXCL5", "IGFL2", "ADAM32", "ADAM18", "IGFL1", "FGF8", "FGF19", "FDF4", "FGF23")




############################################
# Tumour atlas
library(SeuratDisk)
Convert("TICAtlas.h5ad", "TICAtlas.h5seurat")
atlas <- LoadH5Seurat("TICAtlas.h5seurat")

DimPlot(samples.combined, reduction = "umap", label = T)
m0 <- FindConservedMarkers(samples.combined, ident.1 = 0, grouping.var = "orig.ident")
m1 <- FindConservedMarkers(samples.combined, ident.1 = 1, grouping.var = "orig.ident")
m2 <- FindConservedMarkers(samples.combined, ident.1 = 2, grouping.var = "orig.ident")
m3 <- FindConservedMarkers(samples.combined, ident.1 = 3, grouping.var = "orig.ident")
m4 <- FindConservedMarkers(samples.combined, ident.1 = 4, grouping.var = "orig.ident")
m5 <- FindConservedMarkers(samples.combined, ident.1 = 5, grouping.var = "orig.ident")
m6 <- FindConservedMarkers(samples.combined, ident.1 = 6, grouping.var = "orig.ident")
m7 <- FindConservedMarkers(samples.combined, ident.1 = 7, grouping.var = "orig.ident")
m8 <- FindConservedMarkers(samples.combined, ident.1 = 8, grouping.var = "orig.ident")


################3


anchors <- FindTransferAnchors(reference = a, query = samples.combined,
                                        dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = a$cell_type,
                            dims = 1:30)
samples.combined <- AddMetaData(samples.combined, metadata = predictions)


#########

samples.combined <- RenameIdents(samples.combined, `0`="M2-like Macrophages",`1`="M2-like Macrophages",`2`="Fibroblast1",`3`="Malignant",`4`="M1 Macrophages",`5`="M1 Macrophages",`6`="Fibroblast2",
                                 `7`="Malignant",`8`="T / NK cells",`9`="Malignant",`10`="Plasma B cells",`11`="DC",`12`="Proliferative",`13`="Erythrocytes",`14`="pDC")

samples.combined@meta.data[samples.combined@meta.data$patient=="p1",]$condition <- "recurrent"
samples.combined@meta.data[samples.combined@meta.data$patient=="p2",]$condition <- "recurrent"
samples.combined@meta.data[samples.combined@meta.data$patient=="p3",]$condition <- "diagnosis"
samples.combined@meta.data[samples.combined@meta.data$patient=="p4",]$condition <- "initial_treatment"
samples.combined@meta.data[samples.combined@meta.data$sample_ID=="3288",]$condition <- "diagnosis"
samples.combined@meta.data[samples.combined@meta.data$sample_ID=="3288.1",]$condition <- "initial_treatment"
samples.combined@meta.data[samples.combined@meta.data$patient=="p6",]$condition <- "diagnosis"


samples.combined@meta.data[samples.combined@meta.data$patient=="p1",]$sur <- "Upfront_surgery"
samples.combined@meta.data[samples.combined@meta.data$patient=="p2",]$sur <- "Neoadjuvant"
samples.combined@meta.data[samples.combined@meta.data$patient=="p3",]$sur <- "no_surgery"
samples.combined@meta.data[samples.combined@meta.data$patient=="p4",]$sur <- "Neoadjuvant"
samples.combined@meta.data[samples.combined@meta.data$patient=="p5",]$sur <- "Neoadjuvant"
samples.combined@meta.data[samples.combined@meta.data$patient=="p6",]$sur <- "Neoadjuvant"
Idents(samples.combined) <- samples.combined$sub.cluster
samples.combined <- RenameIdents(samples.combined, `0`="M2-like Macrophages",`1`="M2-like Macrophages",`2`="Fibroblast1",`4`="M1 Macrophages",`5`="M1 Macrophages",`6`="Fibroblast2",
                                 `7`="Malignant",`8`="T / NK cells",`9`="Malignant",`10`="Plasma B cells",`11`="DC",`12`="Proliferative",`13`="Plasma B cells",`14`="pDC",
                                 `3_0`="M1 Macrophages", `3_1` = "M2-like Macrophages",`3_2` = "M2-like Macrophages",`3_3` = "Malignant")
