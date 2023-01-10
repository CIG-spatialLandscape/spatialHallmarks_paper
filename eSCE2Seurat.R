##################################################
## Project: Cancer Hallmarks
## Script purpose: Create a Seurat object from an enhanced SingleCellExperiment object
## Date: 22/12/2022
## Author: Sergi Cervilla & Mustafa Sibai
##################################################
library(BayesSpace)
library(Seurat)

#Load spot level single cell experiment object
ST_sce <- readRDS("Desktop/IJC/datasets/Public/P8/RDS/P8_sce.rds") 
#Load sub-spot level single cell experiment object after imputing gene expression
ST_sce.enhanced <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/P8_sce_enhanced.rds") 
#Create Seurat object from enhanced single cell experiment object
STobject_enhanced <- CreateSeuratObject(counts = as.matrix(ST_sce.enhanced@assays@data$SCT), assay = 'SCT', meta.data = as.data.frame(colData(ST_sce.enhanced)))
STobject_enhanced <- SetIdent(STobject_enhanced, value = "spatial.cluster")

## merge the image with STobject_enhanced
STobject.images <- Read10X_Image("Desktop/IJC/datasets/Public/P8/spatial/", image.name = 'tissue_hires_image.png')
STobject.images@assay <- "Spatial"
STobject.images@key <- "P8"
colnames(ST_sce@colData)[3:4] <- c("spot.row","spot.col")

subspot_barcode <- transform(merge(
  x = data.frame(ST_sce@colData)[c("spot", "spot.row", "spot.col")],
  y = cbind(rownames = rownames(data.frame(ST_sce.enhanced@colData)), data.frame(ST_sce.enhanced@colData)),
  by = c("spot.row", "spot.col")), row.names = rownames, rownames = NULL)

subspot_barcode <- tibble::rownames_to_column(subspot_barcode, "subspot_uniq")

colnames(STobject.images@coordinates)[2:3] <- c("spot.row","spot.col")
STobject.images@coordinates <- merge(
  x = subspot_barcode[c("subspot_uniq", "spot.row", "spot.col", "row", "col", "imagerow", "imagecol")],
  y = data.frame(STobject.images@coordinates)[2:3], by = c("spot.row", "spot.col"))

STobject.images@coordinates <- data.frame(STobject.images@coordinates, row.names = "subspot_uniq")
STobject.images@coordinates[1:2] <- NULL

STobject_enhanced@images <- list(P8 = STobject.images)

DefaultAssay(STobject_enhanced) <- "SCT"
STobject_enhanced$spatial.cluster <- as.factor(STobject_enhanced$spatial.cluster)

#Scale data
STobject_enhanced <- ScaleData(STobject_enhanced)
saveRDS("Desktop/IJC/datasets/Public/P8/RDS/P8_eSeurat.rds")