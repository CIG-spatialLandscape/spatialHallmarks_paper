###########################################

#         Olalekan 2021 dataset

##########################################


################################################################################
#                         Load data
################################################################################


#
path_folder <- "Desktop/IJC/datasets/olalekan2021/GSE147082_RAW/"
files <- list.files("Desktop/IJC/datasets/olalekan2021/GSE147082_RAW/")
samples_list <- list()
p <- c(2,5,3,1,4,6)
files <- files[p]
for (i in 1:length(files)) {
  counts <- read.csv(paste0(path_folder, files[i]), row.names = 1)
  
  samples_list[paste0("P",i)] <- Seurat::CreateSeuratObject(counts = counts, project = paste0("P",i,"olalekan"))
}

################################################################################
#                         QC and filtering
################################################################################


samples_list[[1]] <- PercentageFeatureSet(samples_list[[1]], "^MT.", col.name = "percent_mito")
samples_list[[1]] <- PercentageFeatureSet(samples_list[[1]], "^RP[SL]", col.name = "percent_ribo")
VlnPlot(samples_list[[1]], features = c("nFeature_RNA", "nCount_RNA"))

for (i in 1:length(files)) {
  samples_list[[i]] <- PercentageFeatureSet(samples_list[[i]], "^MT.", col.name = "percent_mito")
  samples_list[[i]] <- PercentageFeatureSet(samples_list[[i]], "^RP[SL]", col.name = "percent_ribo")
  VlnPlot(samples_list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo"))
  
  selected_cells <- WhichCells(samples_list[[i]], expression = nFeature_RNA > 600)
  samples_list[[i]] <- subset(samples_list[[i]], cells = selected_cells)
  total_counts_per_cell <- colSums(samples_list[[i]]@assays$RNA@counts)
  mt_genes <- rownames(samples_list[[i]])[grep("^MT.", rownames(samples_list[[i]]))]
  samples_list[[i]]$percent_mito <- colSums(samples_list[[i]]@assays$RNA@counts[mt_genes, ])/total_counts_per_cell
  selected_mt <- WhichCells(samples_list[[i]], expression = percent_mito < 0.1)
  samples_list[[i]] <- subset(samples_list[[i]], cells = selected_mt)
  samples_list[[i]] <- samples_list[[i]][!grepl("^MT.", rownames(samples_list[[i]])),]
  samples_list[[i]] <- samples_list[[i]][!grepl("^RP[SL]", rownames(samples_list[[i]])),]
}

################################################################################
#                         Check batch effects
################################################################################


samples_merged <- merge(samples_list[[1]], y = samples_list[2:6], add.cell.ids = seq(1:length(samples_list)), project = "samples")

samples_merged <- SCTransform(samples_merged)

samples_merged <- RunPCA(samples_merged, verbose = FALSE)
samples_merged <- RunUMAP(samples_merged, reduction = "pca", dims = 1:30)

#clustering
samples_merged <- FindNeighbors(samples_merged, reduction = "pca", dims = 1:30)
samples_merged <- FindClusters(samples_merged, resolution = 0.3)

DimPlot(samples_merged, reduction = "umap", group.by = "orig.ident", repel=TRUE)

saveRDS(samples_merged, "olalekan2021/RDS/merged.rds")

################################################################################
#                         integration
################################################################################


samples_list <- lapply(X = samples_list, FUN = SCTransform) 
#Select variable features among all datasets
features <- SelectIntegrationFeatures(object.list = samples_list, nfeatures = 3000)

samples_list <- PrepSCTIntegration(object.list = samples_list, anchor.features = features)

#select anchors from datasets
samples.anchors <- FindIntegrationAnchors(object.list = samples_list, normalization.method = "SCT",
                                          anchor.features = features)
#rm(patient.list)
#saveRDS(patient.anchors, file="anchors.rds")

#integrate datasets
samples.combined <- IntegrateData(anchorset = samples.anchors, normalization.method = "SCT", features.to.integrate = )
DefaultAssay(samples.combined) <- "integrated"
#dimensional reduction
samples.combined <- RunPCA(samples.combined, verbose = FALSE, seed.use = NULL)
samples.combined <- RunUMAP(samples.combined, reduction = "pca", dims = 1:50)

#clustering
samples.combined <- FindNeighbors(samples.combined, reduction = "pca", dims = 1:50)
samples.combined <- FindClusters(samples.combined, resolution =0.1)

p1 <- DimPlot(samples.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(samples.combined, reduction = "umap", repel = TRUE, label = T)
p1 + p2

saveRDS(samples.combined, "olalekan2021/RDS/integrated.rds")


################################################################################
#                         Manual annotation
################################################################################

markers <- FindAllMarkers(samples.combined, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
conserved_markers1 <- FindConservedMarkers(samples.combined, ident.1 = 1, grouping.var = "orig.ident")

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


CSC <- c("ABCG2", "ALDH1A1", "ALDH1A2", "ALDH1A3", "CD24","CD44", "EPCAM", "KIT", "MYD88", "PROM1","LGR5")
