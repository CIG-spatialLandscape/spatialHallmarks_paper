##################################################
## Project: Cancer Hallmarks
## Script purpose: Create a Seurat object from an enhanced SingleCellExperiment object and compute Hallmark activities
## Author: Sergi Cervilla* & Mustafa Sibai*
##################################################

library(Seurat)
library(dplyr)
library(BayesSpace)


#Load spot level single cell experiment object
ST_sce <- readRDS("") 
#Load sub-spot level single cell experiment object after imputing gene expression
ST_sce.enhanced <- readRDS("") 
#Create Seurat object from enhanced single cell experiment object
STobject_enhanced <- CreateSeuratObject(counts = as.matrix(ST_sce.enhanced@assays@data$SCT), assay = 'SCT', meta.data = as.data.frame(colData(ST_sce.enhanced)))
STobject_enhanced <- SetIdent(STobject_enhanced, value = "spatial.cluster")

## merge the image with STobject_enhanced
STobject.images <- Read10X_Image("", image.name = 'tissue_hires_image.png')
STobject.images@assay <- "Spatial"
STobject.images@key <- ""
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
#`` same as STobject.images@key
STobject_enhanced@images <- list(`` = STobject.images)

DefaultAssay(STobject_enhanced) <- "SCT"
STobject_enhanced$spatial.cluster <- as.factor(STobject_enhanced$spatial.cluster)

#Scale data
STobject_enhanced <- ScaleData(STobject_enhanced)

#table hallmarks and their corresponding genes
H <- read.table("", sep = "\t", header = T)

#transform the data frame to a list
H_features <- list("1"= H[H$H =="H1",1],"2"= H[H$H =="H2",1], "3"= H[H$H =="H3",1],"4"=H[H$H =="H4",1],"5"= H[H$H =="H5",1],"6"= H[H$H =="H6",1],"7"= H[H$H =="H7",1],
                   "8"= H[H$H =="H8",1],"9"= H[H$H =="H9",1],"10"= H[H$H =="H10",1],"11"= H[H$H =="H11",1],"12"= H[H$H =="H12",1], "13"= H[H$H =="H13",1])
#Compute Hallmark activity through AddModuleScore
STobject_enhanced <- AddModuleScore(
  STobject_enhanced,
  H_features,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = NULL,
  name = "H",
  seed = 200,
  search = F
)
#Change names of Hallmarks in the metadata
STobject_enhanced@meta.data[,paste0("H",1:13)] <- scale(STobject_enhanced@meta.data[,paste0("H",1:13)])

#save file (heavy file)
saveRDS(STobject_enhanced, "")

