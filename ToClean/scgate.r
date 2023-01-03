dataset <- read.table("Desktop/IJC/datasets/Single_cell_markers.txt", header = T, sep = "\t")



dataset_filtered<- dataset %>%  filter(str_detect(speciesType,"Human"), str_detect(tissueType, c("Peripheral blood", "Ovary")))



panglao_filtered<- panglao %>%  filter(str_detect(species,"Hs"), `cell type` %in% ct)
ct <- c("B cells", "T cells", "B cells memory", "B cells naive", "Basal cells", "Dendritic cells", "Endothelial cells", "Epithelial cells", "Fibroblasts", "Gamma delta T cells", "Macrophages", "Mast cells", "Monocytes", "NK cells", "Pericytes", "Plasma cells", "T memory cells", "T helper cells")
# convert dataframes to a list of named vectors which is the format for CelliD input
a <- dataset_filtered %>%  
  group_by(`cellName`) %>%  
  summarise(geneset = list(`cellMarker`))

a <- setNames(a$geneset, a$`cellName`)

#remove very short signatures

a <- sapply(a, function(x) {
  j <- c()
  x <- strsplit(x, split = ", ")
  for (i in x) { 
    i <- gsub("-", ".", i)
    j <- append(j, i)
  }
  x <- j
})
a <- a[sapply(a, length) >= 20]

gbm_es <- gsva(AverageExpression(samples.combined)$integrated, a, mx.diff=FALSE)
heatmap(gbm_es)



##################
#scGate
library(scGate)
my_scGate_model <- gating_model(name = "plasmaB", signature = plasmab)
samples.combined <- scGate(data = samples.combined, model = my_scGate_model, assay ="RNA")
my_scGate_model <- gating_model(name = "pDC", signature = pDC)
samples.combined <- scGate(data = samples.combined, model = my_scGate_model, assay ="RNA")
my_scGate_model <- gating_model(name = "NaiveT", signature = naivetcells)
samples.combined <- scGate(data = samples.combined, model = my_scGate_model, assay ="RNA")
my_scGate_model <- gating_model(name = "Fibro", signature = all_gs$Fibroblasts)
samples.combined <- scGate(data = samples.combined, model = my_scGate_model, assay ="RNA")
my_scGate_model <- gating_model(name = "Bcell", signature = bcell)
samples.combined <- scGate(data = samples.combined, model = my_scGate_model, assay ="RNA")
my_scGate_model <- gating_model(name = "Mega", signature = all_gs$Megakaryocytes)
samples.combined <- scGate(data = samples.combined, model = my_scGate_model, assay ="RNA")
my_scGate_model <- gating_model(name = "Tcell", signature = all_gs$`T cells`)
samples.combined <- scGate(data = samples.combined, model = my_scGate_model, assay ="RNA")
my_scGate_model <- gating_model(name = "Tcellcd4", signature = a$`CD4+ cytotoxic T cell`)
samples.combined <- scGate(data = samples.combined, model = my_scGate_model, assay ="RNA")
my_scGate_model <- gating_model(name = "MSCs", signature = MSCs)
samples.combined <- scGate(data = samples.combined, model = my_scGate_model, assay ="RNA")

DimPlot(samples.combined, cols = c(list(Impure = "gray", Pure = "green"))) + theme(aspect.ratio = 1)


my_scGate_model <- gating_model(name = "immune", signature = c("IGKC", "IGHG1", "IGHM"), level = 1)  # initialize model with one positive signature
my_scGate_model <- gating_model(model = my_scGate_model, name = "macrophage", signature = c("MS4A1-",
                                                                                            "FCGR1A-"), level = 2)
my_scGate_model <- gating_model(name = "a", signature = )
samples.combined <- scGate(data = samples.combined, model = my_scGate_model, assay = "integrated")
DimPlot(samples.combined, cols = c(list(Impure = "gray", Pure = "green"))) + theme(aspect.ratio = 1)                                                                      


models.DB <- scGate::get_scGateDB()
my_scGate_model <- models.DB$human$generic$Plasma_cell
names(models.DB$human$generic)

Idents(samples.combined) <- samples.combined$manual1
FeaturePlot(samples.combined, features = c("plasmaB_UCell", "NaiveT_UCell", "pDC_UCell", "Fibro_UCell", "Bcell_UCell", "Mega_UCell", "Tcell_UCell", "Tcellcd4_UCell", "MSCs_UCell")) 


######
my_scGate_model <- gating_model(name = "n", signature = negative)
samples.combined <- scGate(data = samples.combined, model = my_scGate_model, assay ="RNA")
my_scGate_model <- gating_model(name = "p", signature = positive)
samples.combined <- scGate(data = samples.combined, model = my_scGate_model, assay ="RNA")
FeaturePlot(samples.combined, features = c("n_UCell", "p_UCell"))
samples.combined <- FindSubCluster(samples.combined, cluster = "Fibroblast", graph.name = "integrated_snn", resolution = 0.1)
DimPlot(samples.combined, reduction = "umap", group.by = "sub.cluster", label = T)

##################
samples.combined$predicted.celltype
DimPlot(samples.combined, reduction = "umap", group.by = "predicted.celltype", label = T) +
FeaturePlot(samples.combined, features = c("predicted.celltype.score"))


###############
library(SingleR)
hpca.se <- celldex::HumanPrimaryCellAtlasData()

sub_ref <- hpca.se
f <- c()
f <-sort(unique(pred$pruned.labels) ) [table(pred$pruned.labels) < 20]
f <- c(f, sort(unique(pred$pruned.labels) ) [table(pred$pruned.labels) < 20])
sub_ref <- sub_ref[,!sub_ref$label.fine %in% f]


pred <- SingleR(GetAssayData(samples.combined, assay="integrated"), ref = sub_ref, labels = sub_ref$label.fine)
samples.combined$singler <- pred$pruned.labels
DimPlot(samples.combined, reduction = "umap", group.by = "singler", label = T) + NoLegend()
table(pred$pruned.labels)
pred2 <- SingleR(GetAssayData(samples.combined, assay="integrated"), ref = hpca.se, clusters = Idents(samples.combined), labels = hpca.se$label.fine)
samples.combined$singler2 <- pred2$pruned.labels
(pred2)
all.markers <- metadata(pred)$de.genes
a$labels <- pred$pruned.labels
a<-GetAssayData(samples.combined, assay="integrated")
# Beta cell-related markers
library(scater)
plotHeatmap(a, order_columns_by="labels", features=unique(unlist(all.markers$beta)))


##########
sub <- subset(samples.combined, idents = c("Macrophages", "Plasma B cells", "B cells", "T cells", "pDC", "Mast cells")) 

ref.anchors <- FindTransferAnchors(reference = ref, query = sub,
                                        dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = ref.anchors, refdata = ref$cell_type,
                            dims = 1:30)
sub <- AddMetaData(sub, metadata = predictions)
DimPlot(sub, reduction = "umap", group.by = "predicted.id", label = T) + 
FeaturePlot(sub, features = "prediction.score.max")

##########
DimPlot(samples.combined, reduction = "umap", group.by = "all_gs_prediction", label = T)
DimPlot(samples.combined, reduction = "umap", group.by = "ref_cell_gs_prediction", label = T)
DimPlot(samples.combined, reduction = "umap", group.by = "ref_group_gs_prediction", label = T)
DimPlot(samples.combined, reduction = "umap", group.by = "all_gs_prediction2", label = T)
DimPlot(samples.combined, reduction = "umap", group.by = "all_gs_prediction3", label = T)
DimPlot(samples.combined, reduction = "umap", group.by = "scina", label = T)



DimPlot(samples.combined, reduction = "umap", group.by = "all_gs_prediction3", label = T)

f <-sort(unique(samples.combined$all_gs_prediction) ) [table(samples.combined$all_gs_prediction) < 10]
f <- c(f, sort(unique(samples.combined$all_gs_prediction3) ) [table(samples.combined$all_gs_prediction3) < 25])
panglao_all <- panglao %>%  filter(str_detect(species,"Hs"), !str_detect(`cell type`, f))
panglao_all <- panglao %>%  filter(str_detect(species,"Hs"))
panglao_all <- panglao_all[!panglao_all$`cell type` %in% f,]
# convert dataframes to a list of named vectors which is the format for CelliD input
panglao_all <- panglao_all %>%  
  group_by(`cell type`) %>%  
  summarise(geneset = list(`official gene symbol`))
all_gs <- setNames(panglao_all$geneset, panglao_all$`cell type`)

#remove very short signatures
all_gs <- all_gs[sapply(all_gs, length) >= 10]


HGT_all_gs <- RunCellHGT(samples.combined, pathways = all_gs, dims = 1:50, n.features = 200)

# For each cell, assess the signature with the lowest corrected p-value (max -log10 corrected p-value)
all_gs_prediction <- rownames(HGT_all_gs)[apply(HGT_all_gs, 2, which.max)]

# For each cell, evaluate if the lowest p-value is significant
all_gs_prediction_signif <- ifelse(apply(HGT_all_gs, 2, max)>2, yes = all_gs_prediction, "unassigned")

# Save cell type predictions as metadata within the Seurat object
samples.combined$all_gs_prediction3 <- all_gs_prediction_signif

DimPlot(samples.combined, reduction = "umap", group.by = "all_gs_prediction3", label = T)

