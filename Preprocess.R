##################################################
## Project: Cancer Hallmarks
## Script purpose: Preprocess samples, transform to SingleCellExperiment object and run clustering (example of P8)
## Date: 22/12/2022
## Author: Sergi Cervilla & Mustafa Sibai
##################################################

library(STutility)
library(dplyr)
library(Seurat)
library(patchwork)
library(BayesSpace)


#create directory to store figures
if(!file.exists("../P8/figures")) dir.create("../P8/figures") 
#create directory to store RDS files
if(!file.exists("../P8/RDS")) dir.create("../P8/RDS") 


#path of each file
infoTable <- data.frame(samples="../P8/filtered_feature_bc_matrix.h5")
infoTable$spotfiles="../P8/spatial/tissue_positions_list.csv"
infoTable$imgs="../P8/spatial/tissue_hires_image.png"
infoTable$json="../P8/spatial/scalefactors_json.json"

#create Seurat object through STutility, keep genes that have at least 5 counts in the whole tissue
STobject <- InputFromTable(infotable = infoTable, 
                           platform =  "Visium", minUMICountsPerGene = 5)



annotLookup <- read.table("Desktop/IJC/datasets/annotLook.txt", header = T, sep = "\t")

# remove mitochondrial, ribosomal and non-coding genes
protein_genes <- annotLookup$gene_name[annotLookup$gene_type %in% c("protein_coding", "TR_V_gene", "TR_D_gene", "TR_J_gene", "TR_C_gene", "IG_LV_gene", "IG_V_gene", "IG_J_gene", "IG_C_gene" , "IG_D_gene")]

selected_genes <- rownames(STobject)
#keep protein coding genes
selected_genes <- protein_genes[protein_genes %in% selected_genes]
#remove ribosomal genes
selected_genes <- selected_genes[!grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", selected_genes)]
#remove mitochondrial genes
selected_genes <- selected_genes[!grepl("^MT-", selected_genes)]


STobject <- SubsetSTData(STobject, features = selected_genes)


#Load the image in Staffli object to use ManualAnnotation
STobject <- LoadImages(STobject, verbose = TRUE)
ImagePlot(STobject)

#Select spots that are out of the tissue and filter out
STobject<- ManualAnnotation(STobject)
STobject<- SubsetSTData(STobject, labels == "Default")

#Quality Control -> remove spots having a low number of Features
VlnPlot(STobject, features = c("nFeature_RNA", "nCount_RNA"))
ggsave("../P8/figures/VlnPlot_prefiltering.pdf")
STobject<- subset(STobject, nFeature_RNA > 250)
VlnPlot(STobject, features = c("nFeature_RNA", "nCount_RNA"))
ggsave("../P8/figures/VlnPlot_postfiltering.pdf")

#Load image in Seurat object to use Seurat plotting functions
img <- Seurat::Read10X_Image(image.dir = '../P8/spatial/', image.name = "tissue_lowres_image.png")
Seurat::DefaultAssay(object = img) <- 'RNA'
rownames(img@coordinates) <- paste0(rownames(img@coordinates), "_1")
img <- img[colnames(x = STobject)]
STobject[['P8']] <- img
SpatialFeaturePlot(STobject, features = "nFeature_RNA")
ggsave("../P8/figures/nFeature.pdf")

#Normalization of the data returning all genes (SCT)
STobject <- SCTransform(STobject,return.only.var.genes = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.1, assay = "RNA")
saveRDS(STobject, "../P8/RDS/P8.rds")

########## Create SingleCellExperiment object with SCT normalization ###########

spatial_dir <- "../P8/spatial/"

colData <- read.csv(file.path(spatial_dir, "tissue_positions_list.csv"), 
                    header = FALSE)
colnames(colData) <- c("spot", "in_tissue", "row", "col", 
                       "imagerow", "imagecol")
colData$spot <- paste0(colData$spot, "_1")
rownames(colData) <- colData$spot
colData <- colData[colData$in_tissue > 0, ]

#Create SCE object by transfering SCT assay and metadata (Spatial coordinates)
ST_sce <- SingleCellExperiment(assays =list(SCT = as.matrix(GetAssayData(STobject, slot = "data"))))
ST_sce <- SingleCellExperiment(assays =list(SCT = as.matrix(GetAssayData(STobject, slot = "data"))), colData = colData[colnames(ST_sce@assays@data$SCT),])

metadata(ST_sce)$BayesSpace.data <- list()
metadata(ST_sce)$BayesSpace.data$platform <- "Visium"
metadata(ST_sce)$BayesSpace.data$is.enhanced <- FALSE


# Cluster with BayesSpace
set.seed(100)
# Compute top 2000 HVGs and 15 PCA using the SCT data
ST_sce <- spatialPreprocess(ST_sce, platform = "Visium", n.PCs = 15,
                            n.HVGs = 2000, log.normalize = F, assay.type = "SCT")

# Fin q for clusters
ST_sce <- qTune(ST_sce, qs = seq(2, 20), platform = "Visium",d =15)
p <- qPlot(ST_sce)

pdf("../P8/figures/BayesQplot.pdf", width = 8, height = 7)
print(p)
dev.off()

#number of clusters
q = 11
#number of PCA
d = 15

# Run the clustering 
ST_sce <- spatialCluster(ST_sce, q=q, d=d, init.method = "kmeans", platform = "Visium", nrep = 50000, gamma = 3)
ST_sce$spatial.cluster <- as.factor(ST_sce$spatial.cluster)

palette <- RColorBrewer::brewer.pal(q, name = "Paired")
clusterPlot(ST_sce, palette = palette, size=0.05) + labs(title = "Spot-level clustering")
saveRDS(ST_sce, "../P8/RDS/P8_sce.rds")

