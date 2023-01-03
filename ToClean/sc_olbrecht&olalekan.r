
library(Seurat)
library(stringr)


# Load olalekan datasets
path_folder <- "Desktop/IJC/datasets/olalekan2021/GSE147082_RAW/"
files <- list.files("Desktop/IJC/datasets/olalekan2021/GSE147082_RAW/")
samples_list <- list()
p <- c(2,5,3,1,4,6)
files <- files[p]
for (i in 1:length(files)) {
  counts <- read.csv(paste0(path_folder, files[i]), row.names = 1)
  n <- rownames(counts)
  n <- str_replace(n, "HLA\\.", "HLA\\-") 
  n <- str_replace(n, "\\.AS", "\\-AS")
  rownames(counts) <- n 
  samples_list[paste0("P",i)] <- Seurat::CreateSeuratObject(counts = counts, project = paste0("P",i,"olalekan"))
}




#Load oolbercht 
counts <- readRDS("Desktop/IJC/datasets/olbrecht2021/RDS/2095-Olbrecht_counts_matrix.rds")
metadata <- read.csv("Desktop/IJC/datasets/olbrecht2021/RDS/2093-Olbrecht_metadata.csv")

samples_merged <- CreateSeuratObject(counts)
samples_merged$sample_name <- sapply(colnames(counts), function(x) substr(x, start = 18, stop = 25))
samples_merged@meta.data <-left_join(x = samples_merged@meta.data, y = metadata, by="sample_name")
rownames(samples_merged@meta.data) <- colnames(samples_merged)

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

#All data merged
samples_merged2 <- merge(samples_list[[1]], y = samples_list[2:6], add.cell.ids = seq(1:length(samples_list)), project = "samples")
samples_merged_all <- merge(samples_merged, samples_merged2)

samples_merged_all <- NormalizeData(samples_merged_all)
samples_merged_all <- ScaleData(samples_merged_all)
samples_merged_all <- FindVariableFeatures(samples_merged_all)
samples_merged_all <- RunPCA(samples_merged_all, verbose = FALSE)
samples_merged_all <- RunUMAP(samples_merged_all, reduction = "pca", dims = 1:30)
samples_merged_all <- FindNeighbors(samples_merged_all, reduction = "pca", dims = 1:30)
samples_merged_all <- FindClusters(samples_merged_all, resolution = 0.3)

samples_merged_all$orig.ident[samples_merged_all$orig.ident=="SeuratProject"]  <- samples_merged_all$patient_id[samples_merged_all$orig.ident=="SeuratProject"]

DimPlot(samples_merged_all, reduction = "umap", group.by = "orig.ident", repel=TRUE)

#All data integrated
samples_merged_all$orig.ident[samples_merged_all$orig.ident=="SeuratProject"]  <- samples_merged_all$patient_id[samples_merged_all$orig.ident=="SeuratProject"]
samples_list <- SplitObject(samples_merged_all, split.by="orig.ident")

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

saveRDS(samples.combined, "Desktop/IJC/datasets/olbrecht2021/RDS/integrated_olalekan&olbrecht.rds")

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

DimPlot(samples.combined, reduction = "umap", repel=TRUE,label = T) + 
  DimPlot(samples.combined, reduction = "umap", repel=TRUE,label = T, group.by = "orig.ident")


samples.combined <- RenameIdents(samples.combined, `0`="Epithelial", `1`="CAFs", `2`="TAMs", `3`="Epithelial", `4`="CAFs_Immune2", `5`="Tcells",
                                 `6`="Epithelial_Proliferative", `7`="CAFs_Immune1", `8`="PlasmaB", `9`="EC", `10`="Epithelial", `11`="Myofibroblast",
                                 `12`="Bcells", `13`="Fibroblast_Proliferative", `14`="Proliferative_Macrophages", `15`="CAFs", `16`="pDC", `17`="Mast_cells",
                                 `18`="Epithelial", `19`="EC",
                                 `5_2`="Treg", `5_0`="CD4+Tcells", `5_4`="ProliferativeTcells", `5_3`="NK", `5_1`="CD8+Tcells", `2_0`="TAMs", `2_1`="M1 Macrophages",
                                 `2_2`="DC")



############
library(clusterProfiler)
library(pathview)
library(enrichplot)
require(DOSE)


organism = "org.Hs.eg.db"
library(organism, character.only = T)

clust <- 14
gene_list <- m[m$cluster==clust,]$avg_log2FC
names(gene_list) <- m[m$cluster==clust,]$gene
gene_list <- sort(gene_list, decreasing = T)
gse <- gseGO(geneList = gene_list, ont = "BP", pvalueCutoff = 0.05, OrgDb = organism, keyType = "SYMBOL")

dotplot(gse)
g <- pairwise_termsim(gse)
treeplot(g)

clust <- 1
gene_list <- m[m$cluster==clust,]$avg_log2FC
names(gene_list) <- m[m$cluster==clust,]$gene
gene_list <- sort(gene_list, decreasing = T)
gse1 <- gseGO(geneList = gene_list, ont = "BP", pvalueCutoff = 0.05, OrgDb = organism, keyType = "SYMBOL")

dotplot(gse1)

clust <- 4
gene_list <- m[m$cluster==clust,]$avg_log2FC
names(gene_list) <- m[m$cluster==clust,]$gene
gene_list <- sort(gene_list, decreasing = T)
gse4 <- gseGO(geneList = gene_list, ont = "BP", pvalueCutoff = 0.05, OrgDb = organism, keyType = "SYMBOL")

dotplot(gse4)

clust <- 7
gene_list <- m[m$cluster==clust,]$avg_log2FC
names(gene_list) <- m[m$cluster==clust,]$gene
gene_list <- sort(gene_list, decreasing = T)
gse7 <- gseGO(geneList = gene_list, ont = "BP", pvalueCutoff = 0.05, OrgDb = organism, keyType = "SYMBOL")

dotplot(gse7)

clust <- 15
gene_list <- m[m$cluster==clust,]$avg_log2FC
names(gene_list) <- m[m$cluster==clust,]$gene
gene_list <- sort(gene_list, decreasing = T)
gse15 <- gseGO(geneList = gene_list, ont = "BP", pvalueCutoff = 0.05, OrgDb = organism, keyType = "SYMBOL")

dotplot(gse15, )

##############################
clust <- 0
gene_list <- m[m$cluster==clust,]$avg_log2FC
names(gene_list) <- m[m$cluster==clust,]$gene
gene_list <- sort(gene_list, decreasing = T)
gse0 <- gseGO(geneList = gene_list, ont = "BP", pvalueCutoff = 0.05, OrgDb = organism, keyType = "SYMBOL")

dotplot(gse0)
g <- pairwise_termsim(gse0)
treeplot(g)

clust <- 3
gene_list <- m[m$cluster==clust,]$avg_log2FC
names(gene_list) <- m[m$cluster==clust,]$gene
gene_list <- sort(gene_list, decreasing = T)
gse3 <- gseGO(geneList = gene_list, ont = "ALL", pvalueCutoff = 0.05, OrgDb = organism, keyType = "SYMBOL")

dotplot(gse3)
g <- pairwise_termsim(gse3)
treeplot(g)

clust <- 6
gene_list <- m[m$cluster==clust,]$avg_log2FC
names(gene_list) <- m[m$cluster==clust,]$gene
gene_list <- sort(gene_list, decreasing = T)
gse6 <- gseGO(geneList = gene_list, ont = "BP", pvalueCutoff = 0.05, OrgDb = organism, keyType = "SYMBOL")

dotplot(gse6)
g <- pairwise_termsim(gse6)
treeplot(g)

clust <- 10
gene_list <- m[m$cluster==clust,]$avg_log2FC
names(gene_list) <- m[m$cluster==clust,]$gene
gene_list <- sort(gene_list, decreasing = T)
gse10 <- gseGO(geneList = gene_list, ont = "BP", pvalueCutoff = 0.05, OrgDb = organism, keyType = "SYMBOL")

dotplot(gse10)
g <- pairwise_termsim(gse10)
treeplot(g)

clust <- 18
gene_list <- m[m$cluster==clust,]$avg_log2FC
names(gene_list) <- m[m$cluster==clust,]$gene
gene_list <- sort(gene_list, decreasing = T)
gse18 <- gseGO(geneList = gene_list, ont = "BP", pvalueCutoff = 0.05, OrgDb = organism, keyType = "SYMBOL")

dotplot(gse18)
g <- pairwise_termsim(gse18)
treeplot(g)

kegg <- gseKEGG(geneList = gene_list, organism ="human", keyType = "SYMBOL")
head(kegg)


#########################
samples.combined <- FindSubCluster(samples.combined, graph.name = "integrated_snn", cluster = 5)
samples.combined <- FindSubCluster(samples.combined, graph.name = "integrated_snn", cluster = 2, resolution = 0.1)


##########33
m <- FindMarkers(samples.combined, ident.1 = "2_2", only.pos = T)
gene_list <- m$avg_log2FC
names(gene_list) <- rownames(m)
gene_list <- sort(gene_list, decreasing = T)
gse_2 <- gseGO(geneList = gene_list, ont = "BP", pvalueCutoff = 0.05, OrgDb = organism, keyType = "SYMBOL")

dotplot(gse_2)


m <- FindMarkers(samples.combined, ident.1 = "2_0", only.pos = T)
gene_list <- m$avg_log2FC
names(gene_list) <- rownames(m)
gene_list <- sort(gene_list, decreasing = T)
gse_0 <- gseGO(geneList = gene_list, ont = "BP", pvalueCutoff = 0.05, OrgDb = organism, keyType = "SYMBOL")

dotplot(gse_0)

m <- FindMarkers(samples.combined, ident.1 = "2_1", only.pos = T)
gene_list <- m$avg_log2FC
names(gene_list) <- rownames(m)
gene_list <- sort(gene_list, decreasing = T)
gse_1 <- gseGO(geneList = gene_list, ont = "BP", pvalueCutoff = 0.05, OrgDb = organism, keyType = "SYMBOL")

dotplot(gse_1)


m <- FindMarkers(samples.combined, ident.1 = "14", only.pos = T)
gene_list <- m$avg_log2FC
names(gene_list) <- rownames(m)
gene_list <- sort(gene_list, decreasing = T)
gse_14 <- gseGO(geneList = gene_list, ont = "BP", pvalueCutoff = 0.05, OrgDb = organism, keyType = "SYMBOL")

dotplot(gse_14)


samples.combined$test <- "olalekan"

samples.combined$test[colnames(samples.combined) %in% colnames(samples.combined2)] <- samples.combined2$sub.cluster

samples.combined$test[samples.combined$test=="olalekan"] <- samples.combined2$SubTypes
