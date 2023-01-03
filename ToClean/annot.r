################################################################################

#                          Annotation Pipeline

################################################################################

library(Seurat)
library(dplyr)
library(tidyverse)
################################################################################
#                Automatic annotation - Major CellTypes
################################################################################
samples.combined <- readRDS("Desktop/IJC/TFG/sc_integrated.rds")
### Load references databases / datasets 

ref <- readRDS("Desktop/IJC/datasets/OC/RDS/subset.rds") #Paula Nieto - Immune cell atlas subset (ovarian)
#load markers from markers.r

panglao <- read_tsv("Downloads/PanglaoDB_markers_27_Mar_2020.tsv.gz")
#apply filter (strongly recommended)
panglao_filtered <- panglao %>%  filter(str_detect(species,"Hs"), `organ` %in% c("Connective tissue", "Epithelium",
  "Immune system", "Vasculature", "Smooth muscle"), `canonical marker` == 1 ) %>%
  group_by(`cell type`) %>%  
  summarise(geneset = list(`official gene symbol`))
filtered_gs <- setNames(panglao_filtered$geneset, panglao_filtered$`cell type`)
filtered_gs <- filtered_gs[sapply(filtered_gs, length) >= 5]

cellMarker <- read.table("Desktop/IJC/datasets/Summer/Single_cell_markers.txt", header = T, sep = "\t")
#apply filtering (strongly recommended)
cellMarker_filtered <- cellMarker %>%  filter(str_detect(speciesType,"Human"), str_detect(tissueType, c("Peripheral blood", "Ovary")))
cellMarker_filtered <- cellMarker_filtered %>%  
  group_by(`cellName`) %>%  
  summarise(geneset = list(`cellMarker`))
cellMarker_filtered <- setNames(cellMarker_filtered$geneset, cellMarker_filtered$`cellName`)

#remove very short signatures

cellMarker_filtered <- sapply(cellMarker_filtered, function(x) {
  j <- c()
  x <- strsplit(x, split = ", ")
  for (i in x) { 
    i <- gsub("-", ".", i)
    j <- append(j, i)
  }
  x <- j
})




###  GSVA
library(GSVA)
gbm_es_panglao <- gsva(AverageExpression(samples.combined)$integrated, filtered_gs, mx.diff=FALSE,  min.sz=3)
gbm_es_cellMarker <- gsva(AverageExpression(samples.combined)$integrated, cellMarker_filtered, mx.diff=FALSE, min.sz=2)
pheatmap(gbm_es_panglao, angle_col = 315,  rev(paletteer::paletteer_c("grDevices::RdBu", 60) ))
pheatmap(gbm_es_cellMarker, angle_col = 315,  rev(paletteer::paletteer_c("grDevices::RdBu", 60) ))

### SingleR
library(SingleR)
hpca.se <- celldex::HumanPrimaryCellAtlasData()
#immunecellatlas <- celldex::DatabaseImmuneCellExpressionData()

pred_SingleR_clusters_main <- SingleR(GetAssayData(samples.combined, assay="integrated"), ref = hpca.se, clusters = Idents(samples.combined), labels = hpca.se$label.main)
plotScoreHeatmap(pred_SingleR_clusters_main)
pred_SingleR_clusters_fine <- SingleR(GetAssayData(samples.combined, assay="integrated"), ref = hpca.se, clusters = Idents(samples.combined), labels = hpca.se$label.fine)
plotScoreHeatmap(pred_SingleR_clusters_fine)

pred_SingleR_cell_main <- SingleR(GetAssayData(samples.combined, assay="integrated"), ref = hpca.se,  labels = hpca.se$label.main)
samples.combined$singler_main <- pred_SingleR_cell_main$pruned.labels
DimPlot(samples.combined, reduction = "umap", group.by = "singler_main", label = T)

### CelliD
library(CelliD) #

samples.combined <- RunMCA(samples.combined)

HGT_panglao <- RunCellHGT(samples.combined, pathways = filtered_gs, dims = 1:50, n.features = 200)
panglao_gs_prediction <- rownames(HGT_panglao)[apply(HGT_panglao, 2, which.max)]
panglao_gs_prediction_signif <- ifelse(apply(HGT_panglao, 2, max)>2, yes = panglao_gs_prediction, "unassigned")
samples.combined$cellid_panglao <- panglao_gs_prediction_signif
DimPlot(samples.combined, reduction = "umap", group.by = "cellid_panglao", label = T)

#not showing good results
HGT_cellmarker <- RunCellHGT(samples.combined, pathways = cellMarker_filtered, dims = 1:50, n.features = 200)
cellmarker_gs_prediction <- rownames(HGT_panglao)[apply(HGT_cellmarker, 2, which.max)]
cellmarker_gs_prediction_signif <- ifelse(apply(HGT_cellmarker, 2, max)>2, yes = cellmarker_gs_prediction, "unassigned")
samples.combined$cellid_cellmarker <- cellmarker_gs_prediction_signif
DimPlot(samples.combined, reduction = "umap", group.by = "cellid_cellmarker", label = T)


#SCINA
library(SCINA)
results = SCINA(samples.combined@assays$integrated@data, filtered_gs, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.999, rm_overlap=F, allow_unknown=T, log_file='SCINA.log')
samples.combined$scina_panglao <- results$cell_labels
DimPlot(samples.combined, reduction = "umap", group.by = "scina_panglao", label = T)

results = SCINA(samples.combined@assays$integrated@data, cellMarker_filtered, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.999, rm_overlap=F, allow_unknown=T, log_file='SCINA.log')
samples.combined$scina_cellmarker <- results$cell_labels
DimPlot(samples.combined, reduction = "umap", group.by = "scina_cellmarker", label = T)




#Can be further filtered cell types with small amount of cells 

#Take immune cells 
sub <- subset(samples.combined, idents = c("Macrophages", "Plasma B cells", "B cells", "T cells", "pDC", "Mast cells")) #extract immune cells
sub <- subset(samples.combined, idents = c(0,1,3,4,5,8,10,11,12,13,14)) #izat
sub <- subset(samples.combined, idents = c(2,7,9,11,15,16)) #olbrecht
sub <- subset(samples.combined, idents = c(2,5,8,12,14,16, 17)) #integration
### Markers based
#gsva
gbm_es_markers <- gsva(AverageExpression(sub)$integrated, markers, mx.diff=FALSE)
heatmap(gbm_es_markers)
#cellid
HGT_markers <- RunCellHGT(sub, pathways = markers, dims = 1:50, n.features = 200)
markers_gs_prediction <- rownames(HGT_markers)[apply(HGT_markers, 2, which.max)]
markers_gs_prediction_signif <- ifelse(apply(HGT_markers, 2, max)>2, yes = markers_gs_prediction, "unassigned")
sub$cellid_markers <- markers_gs_prediction_signif
DimPlot(sub, reduction = "umap", group.by = "cellid_markers", label = T)
#scina
results = SCINA(sub@assays$integrated@data, markers, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.999,sensitivity_cutoff = 0.1, rm_overlap=F, allow_unknown=T, log_file='SCINA.log')
sub$scina_markers <- results$cell_labels
DimPlot(sub, reduction = "umap", group.by = "scina_markers", label = T)

#Seurat anchors
ref.anchors <- FindTransferAnchors(reference = ref, query = sub,
                                   dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = ref.anchors, refdata = ref$cell_type,
                            dims = 1:30)
sub <- AddMetaData(sub, metadata = predictions)

DimPlot(sub, reduction = "umap", group.by = "predicted.id", label = T) + 
FeaturePlot(sub, features = "prediction.score.max")


## Path enrichment to determine unknown clusters 
samples.combined <- cerebroApp::getMarkerGenes(samples.combined,
                                               groups = c('seurat_clusters'),
                                               assay = "integrated",
                                               organism = "hg"
)

# Get enriched pathways through cerebro
samples.combined <- cerebroApp::getEnrichedPathways(samples.combined,
                                                    databases = c("GO_Biological_Process_2018",
                                                                  "GO_Cellular_Component_2018",
                                                                  "GO_Molecular_Function_2018",
                                                                  "KEGG_2016",
                                                                  "WikiPathways_2016",
                                                                  "Reactome_2016",
                                                                  "Panther_2016",
                                                                  "Human_Gene_Atlas",
                                                                  "Mouse_Gene_Atlas"),
                                                    adj_p_cutoff = 0.01,
                                                    max_terms = 100)

View(samples.combined@misc$enriched_pathways$cerebro_seurat_enrichr$seurat_clusters)


#Name Major subtypes
samples.combined <- RenameIdents(samples.combined, `0`="Epithelial", `1`="Fibroblast", `2`="T cells", `3`="Macrophages", `4`="ESCs (+ Proliferative activity)", `5`="B cell", `6`="Plasma B cells", `7`="Unknown - P6", `8` = "Smooth muscle", `9` = "Endothelial cells", `10` = "pDC", `11`="Mast cells")
samples.combined$MajorCellTypes <- Idents(samples.combined)
DimPlot(samples.combined, reduction = "umap", label = T)

#Subclustering

################## T cells #####################################################
samples.combined <- FindSubCluster(samples.combined, cluster = "T cells", graph.name = "integrated_snn")
DimPlot(samples.combined, reduction = "umap", label = T, group.by = "sub.cluster")
DoHeatmap(samples.combined, )

m <- FindAllMarkers(samples.combined, only.pos = T)
top20 <- m %>% group_by(cluster) %>% top_n(20, avg_log2FC)
DoHeatmap(samples.combined, features = top20$gene)

sub <- subset(samples.combined, idents = "T cells")
sub <- FindSubCluster(sub, cluster = "T cells", graph.name = "integrated_snn")
DimPlot(sub, reduction = "umap", label = T, group.by = "sub.cluster")
m_sub <- FindAllMarkers(sub, only.pos = T)
top20 <- m_sub %>% group_by(cluster) %>% top_n(20, avg_log2FC)
DoHeatmap(sub, features = top20$gene)

FeaturePlot(sub, features = markers$`Naive T cells`)
FeaturePlot(sub, features = markers$`T helper`)
FeaturePlot(sub, features = markers$`Th17 cells`)
FeaturePlot(sub, features = markers$`Proliferative T cells`)
FeaturePlot(sub, features = markers$`CD4 transitional memory `)
DotPlot(sub, features = markers$`Naive memory CD4 T cells`)
FeaturePlot(sub, features = markers$NK)
FeaturePlot(sub, features = markers$`CD8 effector memory `)
FeaturePlot(sub, features = markers$`CD8 pre-exhausted`)
FeaturePlot(sub, features = markers$`CD8 terminally exhausted`)
FeaturePlot(sub, features = markers$`CD8  Cytotoxic T cells`)



my_scGate_model <- gating_model(name = "plasmaB", signature = plasmab)

gbm_es_panglao <- gsva(AverageExpression(sub)$integrated, cellMarker_filtered, mx.diff=FALSE)
heatmap(gbm_es_panglao)

Idents(samples.combined)
samples.combined <- RenameIdents(samples.combined, `T cells_0` = "CD4+ Naive T cells", `T cells_2` = "T helper", `T cells_1` = "CD8+ T cells", `T cells_3` = "CD8+ T cells - NK")
paste0(top20[top20$cluster=="T cells_1",]$gene[1:10], collapse = ", ")

################## Macrophages #####################################################
samples.combined <- FindSubCluster(samples.combined, cluster = "Macrophages", graph.name = "integrated_snn")
top20 <- m %>% group_by(cluster) %>% top_n(20, avg_log2FC)
DoHeatmap(samples.combined, features = top20$gene)

sub <- subset(samples.combined, idents = c("Macrophages_0", "Macrophages_3", "Macrophages_1", "Macrophages_2"))
DimPlot(sub, reduction = "umap", label = T)
sub <- RenameIdents(sub, `Macrophages_0` = "Macrophages_1")
DotPlot(sub, features = markers$`Macrophages spp1`)
DotPlot(sub, features = markers$`Macro. and mono. prolif`)
DotPlot(sub, features = markers$`TAMs C1QC`)
DotPlot(sub, features = markers$Monocytes)
samples.combined <- RenameIdents(samples.combined, `Macrophages_0` = "TAMs", `Macrophages_1`="TAMs", `Macrophages_2` = "NR1H2+ Macrophages", `Macrophages_3`="M1 Macrophages")

################## Fibroblast #####################################################


samples.combined <- FindSubCluster(samples.combined, cluster = "Fibroblast", graph.name = "integrated_snn", resolution = 0.2)
DimPlot(samples.combined, reduction = "umap", label = T, group.by = "sub.cluster")
Idents(samples.combined) <- samples.combined$sub.cluster
m <- FindAllMarkers(samples.combined, only.pos = T)
top20 <- m %>% group_by(cluster) %>% top_n(20, avg_log2FC)
DoHeatmap(samples.combined, features = top20$gene)


sub <- subset(samples.combined, ident = "Fibroblast")
Idents(sub) <- sub$sub.cluster
m <- FindAllMarkers(sub, only.pos = T)
top20 <- m %>% group_by(cluster) %>% top_n(20, avg_log2FC)
DoHeatmap(sub, features = top20$gene)
