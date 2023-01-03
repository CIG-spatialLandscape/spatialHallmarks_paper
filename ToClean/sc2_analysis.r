
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)


setwd('Desktop/IJC/datasets')

counts.sample = read.table("data.raw.matrix.txt", header = T, row.names = 1)
sc <- Seurat::CreateSeuratObject(counts = counts.sample)
patients <- stri_sub(colnames(sc), -2, -1)
sc$patients <- patients
VlnPlot(sc, features = "nCount_RNA")
VlnPlot(sc, features = "nFeature_RNA")



total_counts_per_cell <- colSums(sc@assays$RNA@counts)
mt_genes <- rownames(sc)[grep("^MT-", rownames(sc))]
sc$percent_mito <- colSums(sc@assays$RNA@counts[mt_genes, ])/total_counts_per_cell
selected_mt <- WhichCells(sc, expression = percent_mito < 0.1)
sc <- subset(sc, cells = selected_mt)
sc <- sc[!grepl("^MT-", rownames(sc)),]
sc <- sc[!grepl("^RP[SL]", rownames(sc)),]

sc <- SCTransform(sc)

sc <- RunPCA(sc, verbose = FALSE)
sc <- RunUMAP(sc, reduction = "pca", dims = 1:30)

#clustering
sc <- FindNeighbors(sc, reduction = "pca", dims = 1:30)
sc <- FindClusters(sc, resolution = 0.3)

DimPlot(sc, reduction = "umap", group.by = "patients", repel=TRUE)

Idents(prostate_sc) <- prostate_sc$cellid_pred
samples_list <- SplitObject(sc, split.by = "patients")

features <- SelectIntegrationFeatures(samples_list , nfeatures = 3000)

samples_list <- PrepSCTIntegration(object.list = samples_list, anchor.features = features)

#select anchors from datasets
samples.anchors <- FindIntegrationAnchors(object.list = samples_list, normalization.method = "SCT",
                                          anchor.features = features)
#rm(patient.list)
#saveRDS(patient.anchors, file="anchors.rds")

#integrate datasets
samples.combined <- IntegrateData(anchorset = samples.anchors, normalization.method = "SCT")

samples.combined$cellid_pred[samples.combined$cellid_pred=="B cells memory"]<-"B cells"
samples.combined$cellid_pred[samples.combined$cellid_pred=="Basophils"]<-"Chondrocytes"

prostate_sc$cellid_pred2 <- prostate_sc$cellid_pred
prostate_sc$cellid_pred2[prostate_sc$cellid_pred2=="unassigned"]<-prostate_sc$seurat_clusters[prostate_sc$cellid_pred2=="unassigned"]
sc$cellid_pred[sc$cellid_pred=="Basophils"]<-"Basophils"