library(data.table)
data_to_write_out <- as.data.frame(as.matrix(OC.st@assays$Spatial@counts))
fwrite(x = data_to_write_out, row.names = T,file = "matrix_public.csv")
data_to_write_out <- as.data.frame(as.matrix(OC.st@images$OC@coordinates))
fwrite(x = data_to_write_out, row.names = T,file = "matrix_coordinates_public.csv")

d <- OC.st@images$OV4A@coordinates[,4:5]
distance_m <- dist(d)

prediction <- as.data.frame(as.matrix(OC.st@assays$prediction@data))

prediction[,3]
                                      
value <- 0
for (spot in colnames(distances)) {
   
    sum(prediction[3, ]/distances[spot,])
    
}

library(SeuratDisk)
SaveH5Seurat(samples.combined, "Desktop/IJC/datasets/olbrecht2021/RDS/sc_tumors.h5Seurat", overwrite = T)
Convert("Desktop/IJC/datasets/olbrecht2021/RDS/sc_tumors.h5Seurat", dest= "h5ad", overwrite = T)

SaveH5Seurat(samples.combined, "Desktop/IJC/datasets/olbrecht2021/RDS/sc_integrated.h5Seurat", overwrite = T)
Convert("Desktop/IJC/datasets/olbrecht2021/RDS/sc_integrated.h5Seurat", dest= "h5ad", overwrite = T)

m <- m/rowSums(m)
#m <- m/apply(m,2,max)

m1 <- sapply(m, function(x) {
  x / max(x)
})
apply(m,1,max)




##################3
library(dplyr, BayesSpace)
library(Seurat)

sce <- readRDS("Desktop/OV_4A_enhanced.rds")

OC.st <- as.Seurat(sce, counts = "SCT", data = "SCT")


### Try to label transfer on the enhanced object
names(sce@assays) <- "logcounts"

OC.st <- Seurat::CreateSeuratObject(counts=logcounts(sce),
                                               assay='logcounts',
                                               meta.data=as.data.frame(colData(sce)))

OC.st <- Seurat::SetIdent(OC.st, value = "spatial.cluster")

OC.st <- RenameAssays(OC.st, `originalexp`="SCT")


## Scale data
OC.st@assays$SCT@scale.data <-
  OC.st@assays$SCT@data %>% as.matrix %>% t %>% scale %>% t

## merge the image with OV.enhanced.st
NASH.images <- Read10X_Image("/home/malsibai/rnaseq/data/ST/HCC/processed/HCC_Output/HCC3M265/spatial/")
NASH.images@assay <- "Spatial"
NASH.images@key <- "NASH_"
image <- OC@images$OV4A
colnames(image@coordinates)[2:3] <- c("spot.row","spot.col")

subspot_barcode <- transform(merge(
  x = data.frame(image@coordinates)[,c("spot.row", "spot.col")],
  y = cbind(rownames = rownames(data.frame(sce@colData)), data.frame(sce@colData)),
  by = c("spot.row", "spot.col")), row.names = rownames, rownames = NULL)

subspot_barcode <- tibble::rownames_to_column(subspot_barcode, "subspot_uniq")

colnames(NASH.images@coordinates)[2:3] <- c("spot.row","spot.col")
image@coordinates <- merge(
  x = subspot_barcode[c("subspot_uniq", "spot.row", "spot.col", "row", "col", "imagerow", "imagecol")],
  y = data.frame(image@coordinates)[2:3], by = c("spot.row", "spot.col"))

image@coordinates <- data.frame(image@coordinates, row.names = "subspot_uniq")
image@coordinates[1:2] <- NULL

# now input the images into NASH.Enhanced.st seurat object

NASH.enhanced.st@images <- list(NASH = NASH.images)
OC.st@images <- list(OC=image)

DefaultAssay(OC.st) <- "SCT"
SpatialDimPlot(OC.st, cols = palette.mod , pt.size.factor = 0.4)

NASH.enhanced.st$spatial.cluster <- as.factor(NASH.enhanced.st$spatial.cluster)

palette.mod <- c("#E41A1C", "#377EB8", "#4DAF4A", "#999999", "#984EA3", "#FF7F00", "#A65628", "#FFFF33", "#F781BF")
p <- SpatialDimPlot(NASH.enhanced.st, cols = palette.mod , pt.size.factor = 0.6)

pdf(file = "./EDA_NASH/HCC_Enhanced_BayesSpace_W_HEBack.pdf")
print(p)
dev.off()



##################3

DefaultAssay(OC.st) <- "SCT"
DefaultAssay(samples.combined) <- "SCT"

anchors <- FindTransferAnchors(reference = samples.combined, query = OC.st, normalization.method = "SCT", n.trees=1000, reduction = "cca")

predictions.assay <- TransferData(anchorset = anchors, refdata = Idents(samples.combined), prediction.assay = TRUE,
                                  weight.reduction = "cca", dims = 1:30, n.trees = 1000)

OC.st[["prediction"]] <- predictions.assay
DefaultAssay(OC.st) <- "prediction"
p <- SpatialFeaturePlot(object = OC.st, features = rownames(OC.st), , ncol = 3, pt.size.factor = 0.6)
pdf(file = "Desktop/IJC/datasets/olbrecht2021/figures/4A_proli.pdf", width = 30, height = 30)
plot(p)
dev.off()
