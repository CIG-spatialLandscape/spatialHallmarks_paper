
sc_6 <- readRDS("integrated_seurat.rds")
sc_13 <- readRDS("Single-Cell/sc_combined.rds")

sc_6[["mca"]] <- NULL
sc_13[["mca"]] <- NULL

all.list <- c(sc_6, sc_13)
features <- SelectIntegrationFeatures(object.list = all.list, nfeatures = 3000)

all.list <- PrepSCTIntegration(object.list = all.list, anchor.features = features)

#select anchors from datasets
all.anchors <- FindIntegrationAnchors(object.list = all.list, normalization.method = "SCT",
                                          anchor.features = features)
#rm(patient.list)
#saveRDS(patient.anchors, file="anchors.rds")

#integrate datasets
all.combined <- IntegrateData(anchorset = all.anchors, normalization.method = "SCT")
DefaultAssay(all.combined) <- "integrated"
#dimensional reduction
all.combined <- RunPCA(all.combined, verbose = FALSE)
all.combined <- RunUMAP(all.combined, reduction = "pca", dims = 1:30)

#clustering
all.combined <- FindNeighbors(all.combined, reduction = "pca", dims = 1:30)
all.combined <- FindClusters(all.combined, resolution =0.3)

p1 <- DimPlot(all.combined, reduction = "umap", group.by = "sub.cluster", label=T)
p2 <- DimPlot(all.combined, reduction = "umap", repel = TRUE)
p1 + p2
Idents(all.combined) <- all.combined$cellid_pred
all.combined <- FindSubCluster(all.combined, "unassigned", "integrated_snn", resolution = 0.1)

all.combined$sub.cluster[all.combined$sub.cluster=="Natural killer T cells"]<-"NK cells"
all.combined$sub.cluster[all.combined$sub.cluster=="Megakaryocytes"]<-"Mast cells"
all.combined$sub.cluster[all.combined$sub.cluster=="Adipocytes"]<-"unassigned_2"
all.combined$sub.cluster[all.combined$sub.cluster=="T helper cells"]<-"T cells"
all.combined$sub.cluster[all.combined$sub.cluster=="B cells naive"]<-"B cells"






j = 7
for (i in 1:length(samples_list1)) {
  samples_list[[j]] <- samples_list1[[i]]
  j = j+1
  print(j)
}
