
#author: Sergi Cervulla
#date: 17/2/2022

#### Packages ####
library(STutility)



### Load data
#create table with the file paths (STutlity format)
infoTable <- data.frame(samples="Desktop/IJC/datasets/IGTP/4A/filtered_feature_bc_matrix.h5",
                        spotfiles="Desktop/IJC/datasets/IGTP/4A/images/tissue_positions_list.csv",
                        imgs="Desktop/IJC/datasets/IGTP/4A/images/tissue_hires_image.png",
                        json="Desktop/IJC/datasets/IGTP/4A/images/scalefactors_json.json")

#create Seurat object with filtering
OC_4A <- InputFromTable(infotable = infoTable, 
                        min.gene.count = 100, 
                        min.gene.spots = 5,
                        min.spot.count = 500,
                        platform =  "Visium")

‘imager’, ‘raster’, ‘magick’, ‘imagerExtra’, ‘ggiraph’, ‘spdep’ are not available for package ‘STutility’



#count matrix
expr.data <- Seurat::Read10X_h5(filename = "Desktop/IJC/datasets/IGTP/4A/filtered_feature_bc_matrix.h5")


OC.st <- Seurat::CreateSeuratObject(counts = expr.data, project = 'oc', assay = 'Spatial', na.rm=T)

img <- Seurat::Read10X_Image(image.dir = 'Desktop/IJC/datasets/IGTP/4A/images/')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = OC.st)]
OC.st[['OC']] <- img


### QC
p1 <- ggplot() +
  geom_histogram(data = OC.st[[]], aes(nFeature_Spatial), fill = "red", alpha = 0.7, bins = 50) +
  ggtitle("Unique genes per spot")

p2 <- ggplot() +
  geom_histogram(data = OC.st[[]], aes(nCount_Spatial), fill = "red", alpha = 0.7, bins = 50) +
  ggtitle("Total counts per spots")

gene_attr <- data.frame(nUMI = Matrix::rowSums(OC.st@assays$Spatial@counts), 
                        nSpots = Matrix::rowSums(OC.st@assays$Spatial@counts > 0))
p3 <- ggplot() +
  geom_histogram(data = gene_attr, aes(nUMI), fill = "red", alpha = 0.7, bins = 50) +
  scale_x_log10() +
  ggtitle("Total counts per gene (log10 scale)")

p4 <- ggplot() +
  geom_histogram(data = gene_attr, aes(nSpots), fill = "red", alpha = 0.7,  bins = 50) +
  ggtitle("Total spots per gene")

p <- (p1 - p2)/(p3 - p4)


### Filtering
# remove genes with fewer than 10 counts across the whole dataset 
# remove mitochondrial protein coding genes as well as all non-coding RNA biotypes (“antisense,” “lincRNA,” “pseudogenes”)
gene_list <- read.table("Desktop/annotLookup.xls", header = T)
protein_genes <- gene_list$external_gene_name[gene_list$gene_biotype=="protein_coding"]
selected_genes <- rownames(OC.st)[rowSums(OC.st) > 5]
selected_genes <- protein_genes[protein_genes %in% selected_genes]
selected_genes <- selected_genes[!grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", selected_genes)]
selected_genes <- selected_genes[!grepl("^MTRNR", selected_genes)]

OC.st <- OC.st[selected_genes,]  
mito.genes <- grep("^MT-", rownames(expr.data), value = TRUE)


ggplot(a, aes(x=b, y=log(val))) + geom_violin() + geom_jitter(aes(size=0.1)) + scale_size_area(max_size =0.1)


#### SCTtransform
OC.st <- SCTransform(OC.st,return.only.var.genes = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.1, assay = "Spatial")

OC.st <- RunPCA(OC.st, assay = "SCT", verbose = FALSE)

pca=OC.st@reductions$pca
eival= (pca@stdev)^2
varexplaiend = eival/sum(eival)
plot(1:50, cumsum(varexplaiend), type="l")

OC.st <- FindNeighbors(OC.st, reduction = "pca", dims = 1:50)
OC.st <- FindClusters(OC.st, verbose = FALSE, resolution = 0.7)
OC.st <- RunUMAP(OC.st, reduction = "pca", dims = 1:50)

SpatialDimPlot(OC.st, repel=TRUE,label = T)

OC.st_3 <- RunNMF(OC.st, nfactors = 15)

OC.st_3 <- AddMetaData(OC.st_3, as.data.frame(OC.st_3@reductions$NMF@cell.embeddings))

p <- SpatialFeaturePlot(OC.st_3, features = colnames(OC.st_3@meta.data)[6:26], ncol=3)
pdf("Desktop/IJC/datasets/IGTP/4A/figures/factors15_norb.pdf", width = 10, height = 30)
plot(p)
dev.off()


FactorGeneLoadingPlot(OC.st_3, factor = "7", topn = 30)
