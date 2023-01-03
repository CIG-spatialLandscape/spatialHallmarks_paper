################################
# Author: Sergi Cervilla 
# Reproducible analysis of single cell data
################################


library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)


setwd('Desktop/IJC/datasets')




#samples_paths = c("Single-Cell/GSM4089151_P1_gene_cell_exprs_table.txt", "Single-Cell/GSM4089152_P2_gene_cell_exprs_table.txt",
                  #"Single-Cell/GSM4089153_P3_gene_cell_exprs_table.txt", "", "", "") #vector containing the file for each patient

samples_list = readRDS("Single-Cell/patient.list.rds")

#Filtering
for (i in 1:length(samples_list)) {
  #counts.sample = read.table(samples_paths[i], header = T, row.names = 1)
  #counts.sample <- counts.sample[, -1]
  #rownames(x = counts.sample ) <- make.unique(names = sapply(X = strsplit(x = rownames(x = counts.sample), split = '__'), FUN = '[', 1))
  #project_name = paste0('sample#', i)
  #samples_list[i] <- Seurat::CreateSeuratObject(counts = counts.patient1, project = project_name)
  
  selected_cells <- WhichCells(samples_list[[i]], expression = 7000 > nFeature_RNA & nFeature_RNA > 500)
  samples_list[[i]] <- subset(samples_list[[i]], cells = selected_cells)
  
  total_counts_per_cell <- colSums(samples_list[[i]]@assays$RNA@counts)
  mt_genes <- rownames(samples_list[[i]])[grep("^MT-", rownames(samples_list[[i]]))]
  samples_list[[i]]$percent_mito <- colSums(samples_list[[i]]@assays$RNA@counts[mt_genes, ])/total_counts_per_cell
  selected_mt <- WhichCells(samples_list[[i]], expression = percent_mito < 0.1)
  samples_list[[i]] <- subset(samples_list[[i]], cells = selected_mt)
  samples_list[[i]] <- samples_list[[i]][!grepl("^MT-", rownames(samples_list[[i]])),]
  samples_list[[i]] <- samples_list[[i]][!grepl("^RP[SL]", rownames(samples_list[[i]])),]
}



################# Check batch effect ########################

samples_merged <- merge(samples_list[[1]], y = samples_list[2:6], add.cell.ids = seq(1:length(samples_list)), project = "samples")

samples_merged <- SCTransform(samples_merged)

samples_merged <- RunPCA(samples_merged, verbose = FALSE)
samples_merged <- RunUMAP(samples_merged, reduction = "pca", dims = 1:30)

#clustering
samples_merged <- FindNeighbors(samples_merged, reduction = "pca", dims = 1:30)
samples_merged <- FindClusters(samples_merged, resolution = 0.3)

DimPlot(samples_merged, reduction = "umap", group.by = "orig.ident", repel=TRUE)

#############################################################

###################### Integration ##########################

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
samples.combined <- IntegrateData(anchorset = samples.anchors, normalization.method = "SCT")
DefaultAssay(samples.combined) <- "integrated"
#dimensional reduction
samples.combined <- RunPCA(samples.combined, verbose = FALSE)
samples.combined <- RunUMAP(samples.combined, reduction = "pca", dims = 1:30)

#clustering
samples.combined <- FindNeighbors(samples.combined, reduction = "pca", dims = 1:30)
samples.combined <- FindClusters(samples.combined, resolution =0.3)

p1 <- DimPlot(patient.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(samples.combined, reduction = "umap", repel = TRUE)
p1 + p2

saveRDS(samples.combined, "samples.combined.rds")

#########################Harmony############################

library(harmony)
samples_merged <- readRDS("merged.rds")
samples_merged <- samples_merged %>% RunHarmony("orig.ident")

#############################################################

################### Manual annotation #######################

markers <- FindAllMarkers(samples.combined, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
markers <- markers[order(markers$avg_log2FC, decreasing = T), ]
top20 <- patient.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)

FeaturePlot(samples.combined, features = , label=T)
VlnPlot(samples.combined, features = )
DoHeatmap(samples.combined, features = top10, label = TRUE)


pie(table(panglao[panglao$`official gene symbol` %in% markers[markers$cluster==14,]$gene,]$organ))
pie(tail(sort(table(panglao[panglao$`official gene symbol` %in% markers[markers$cluster==14,]$gene,]$`cell type`)), n=10))

#############################################################

################### Automatic annotation ####################
library(CelliD)
library(stringr)
library(readr)

panglao <- read_tsv("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz")

features <- markerss$gene
features <- unique(features)
samples.combined <- RunMCA(samples.combined)

panglao_filtered<- panglao %>%  filter(str_detect(species,"Hs"), `organ` %in% c("Connective tissue", "Epithelium",
                                                                             "Immune system", "Vasculature"))

panglao_filtered<- panglao %>%  filter(str_detect(species,"Hs"), `cell type` %in% ct)
ct <- c("B cells", "T cells", "B cells memory", "B cells naive", "Basal cells", "Dendritic cells", "Endothelial cells", "Epithelial cells", "Fibroblasts", "Gamma delta T cells", "Macrophages", "Mast cells", "Monocytes", "NK cells", "Pericytes", "Plasma cells", "T memory cells", "T helper cells")
# convert dataframes to a list of named vectors which is the format for CelliD input
panglao_filtered <- panglao_filtered %>%  
  group_by(`cell type`) %>%  
  summarise(geneset = list(`official gene symbol`))
all_gs <- setNames(panglao_filtered$geneset, panglao_filtered$`cell type`)

#remove very short signatures
all_gs <- all_gs[sapply(all_gs, length) >= 10]

HGT_all_gs <- RunCellHGT(samples.combined, pathways = all_gs, dims = 1:20)

all_gs_prediction <- rownames(HGT_all_gs)[apply(HGT_all_gs, 2, which.max)]

all_gs_prediction_signif <- ifelse(apply(HGT_all_gs, 2, max)>2, yes = all_gs_prediction, "unassigned")


samples.combined$cellid_pred <- all_gs_prediction_signif
#############################################################
unique(samples.combined$cellid_pred)
Idents(samples.combined) <- samples.combined$seurat_clusters
p1 <- DimPlot(samples.combined, reduction = "umap", group.by = "cellid_pred", label=TRUE)
p2 <- DimPlot(samples.combined, reduction = "umap", group.by = "patients", label=TRUE)
wrap_plots(p1,p2)
pdf(file="plot_merged.pdf", width=15, height=10)
plot(p)
dev.off()


samples.combined <- IntegrateData(anchorset = samples.anchors, normalization.method = "SCT")
DefaultAssay(patient.combined) <- "integrated"
#dimensional reduction
sc <- RunPCA(sc, verbose = FALSE)
sc <- RunUMAP(sc, reduction = "pca", dims = 1:30)

#clustering
sc <- FindNeighbors(sc, reduction = "pca", dims = 1:30)
sc <- FindClusters(sc, resolution =0.5)

p1 <- DimPlot(patient.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(sc, reduction = "umap", repel = TRUE)
p1 + p2



##############################################################
sc$cellid_pred[sc$cellid_pred=="unassigned"]<-sc$seurat_clusters[sc$cellid_pred=="unassigned"]
sc$cellid_pred[sc$cellid_pred=="Basophils"]<-"Basophils"
sc$cellid_pred[sc$cellid_pred=="Basophils"]<-"Mast cells"
sc$cellid_pred[sc$cellid_pred=="Red pulp macrophages"]<-"Macrophages"
sc$cellid_pred[sc$cellid_pred=="Nuocytes"]<-"T cells"
sc$cellid_pred[sc$cellid_pred=="Endothelial cells (blood brain barrier)"]<-"Endothelial cells"
sc$cellid_pred[sc$cellid_pred=="Adipocyte progenitor cells"]<-"Pericytes"
sc$cellid_pred[sc$cellid_pred=="Mesothelial cells"]<-"4"
sc$cellid_pred[sc$cellid_pred=="1"]<-"Unknown 1"
sc$cellid_pred[sc$cellid_pred=="10"]<-"Unknown 10"
sc$cellid_pred[sc$cellid_pred=="11"]<-"Unknown 11"
sc$cellid_pred[sc$cellid_pred=="12"]<-"Unknown 12"
sc$cellid_pred[sc$cellid_pred=="13"]<-"Unknown 13"
sc$cellid_pred[sc$cellid_pred=="14"]<-"Endothelial cells"
sc$cellid_pred[sc$cellid_pred=="15"]<-"Unknown 15"
sc$cellid_pred[sc$cellid_pred=="21"]<-"Fibroblasts"
sc$cellid_pred[sc$cellid_pred=="18"]<-"B cells"
sc$cellid_pred[sc$cellid_pred=="3"]<-"T cells"
sc$cellid_pred[sc$cellid_pred=="4"]<-"Unknown 4"
sc$cellid_pred[sc$cellid_pred=="5"]<-"Unknown 5"
sc$cellid_pred[sc$cellid_pred=="6"]<-"Unknown 6"
sc$cellid_pred[sc$cellid_pred=="7"]<-"Unknown 7"
sc$cellid_pred[sc$cellid_pred=="8"]<-"Unknown 8"
sc$cellid_pred[sc$cellid_pred=="9"]<-"Unknown 9"
sc$cellid_pred[sc$cellid_pred=="2"]<-"Unknown 2"
sc$cellid_pred[sc$cellid_pred=="16"]<-"Unknown 16"
sc$cellid_pred[sc$cellid_pred=="17"]<-"Unknown 17"
sc$cellid_pred[sc$cellid_pred=="19"]<-"Unknown 19"
sc$cellid_pred[sc$cellid_pred=="20"]<-"Endothelial cells"
sc$cellid_pred[sc$cellid_pred=="22"]<-"Unknown 22"
sc$cellid_pred[sc$cellid_pred=="23"]<-"Endothelial cells"
sc$cellid_pred[sc$cellid_pred=="bbbbbbbbbbb"]<-"Endothelial cells"

sc$cellid_pred[sc$cellid_pred=="Unknown 22"]<-"Unknown 5"
sc$cellid_pred[sc$cellid_pred=="Unknown 19"]<-"Unknown 5"
#sc$cellid_pred[sc$cellid_pred=="Unknown 8"]<-"Endothelial cells"
sc$cellid_pred[sc$cellid_pred=="Unknown 9"]<-"Endothelial cells"
m <- FindMarkers(sc, ident.1 = 2)
