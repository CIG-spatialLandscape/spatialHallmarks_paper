##################################################
## Project: Cancer Hallmarks
## Script purpose: Create files to plot high resolution images with enhanced BayesSpace spots
## Date: 22/12/2022
## Author: Sergi Cervilla* & Mustafa Sibai*
##################################################

## Auxiliary functions
#sort clusters from Cancer(1) to TME(5)
new_labels <- function(original, order){
  labels <- original
  #Cancer (lowest)
  labels[original==order[1]] <- "1"
  labels[original==order[2]] <- "2"
  labels[original==order[3]] <- "3"
  labels[original==order[4]] <- "4"
  #TME (highest)
  labels[original==order[5]] <- "5"
  return(labels)
}


source("../utils/PlottingMod.R")

sample <- "M3"
#load single cell experiment object at spot level
ST_sce <- readRDS("")
#load single cell experiment object at sub-spot level
ST_sce.enhanced <- readRDS("")

## merge the image with STobject.enhanced
#select high resolution image (spatial folder)
STobject.images <- Read10X_Image("", image.name = 'tissue_hires_image.png')
STobject.images@assay <- "Spatial"
STobject.images@key <- sample
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

mat <- matrix(0, nrow=1, ncol=ncol(ST_sce.enhanced))
colnames(mat) <- colnames(ST_sce.enhanced)
rownames(mat) <- "dummy"
STobject <- CreateSeuratObject(mat)

STobject@images <- list(obj = STobject.images)

STobject@images[[1]]@scale.factors$lowres = STobject@images[[1]]@scale.factors$hires
STobject$test <- NA


################################################################################
hires <- SpatialDimPlot(STobject, pt.size.factor = 2, image.alpha = 1) 

#load output from estimate
estimate <- read.table("", sep = "\t", skip = 2, header = T)
estimate <- t(estimate[3, 3:ncol(estimate)])

ST_sce.enhanced <- ST_sce.enhanced[,-1]
labels <- kmeans(estimate, centers = 5)$cluster

labels_order <- order(sapply(1:5, function(x){
  mean(estimate[labels==x])
}))

ST_sce.enhanced$clusters <- new_labels(labels,labels_order)

v <- .make_triangle_subspots(colData(ST_sce.enhanced), fill = "clusters")

if (cor(ST_sce$imagerow, ST_sce$spot.row) > 0.99 ) {
  message("no re-oriented image")
  v$imagerow <- sapply(v$x.vertex, function(x) {
    m <- (x - min(v$x.vertex))/(max(v$x.vertex)-min(v$x.vertex)) * (max(hires$data$imagecol)-min(hires$data$imagecol)) + min(hires$data$imagecol)
  })
  
  v$imagecol <- sapply(v$y.vertex, function(x) {
    m <- (x - min(v$y.vertex))/(max(v$y.vertex)-min(v$y.vertex)) * (max(hires$data$imagerow)-min(hires$data$imagerow)) + min(hires$data$imagerow)
  })
  
} else {
  message("re-oriented image")
  v$x.vertex <- sapply(-v$x.vertex, function(x) {
    m <- (x - min(-v$x.vertex))/(max(-v$x.vertex)-min(-v$x.vertex)) * (max(v$x.vertex)-min(v$x.vertex)) + min(v$x.vertex)
  })
  v$imagecol <- sapply(v$x.vertex, function(x) {
    m <- (x - min(v$x.vertex))/(max(v$x.vertex)-min(v$x.vertex)) * (max(hires$data$imagecol)-min(hires$data$imagecol)) + min(hires$data$imagecol)
  })
  v$y.vertex <- abs(v$y.vertex)
  v$imagerow <- sapply(v$y.vertex, function(x) {
    m <- (x - min(v$y.vertex))/(max(v$y.vertex)-min(v$y.vertex)) * (max(hires$data$imagerow)-min(hires$data$imagerow)) + min(hires$data$imagerow)
  })
}

#apply image shifts to align with H&E image
xshift <- 20
yshift <- 100
hires + geom_polygon(data=v, 
                     aes_(x=~imagerow+xshift, y=~imagecol+yshift, group=~spot, fill=~factor(fill))) 

#apply image shifts
v$imagerow <- v$imagerow + xshift
v$imagecol <- v$imagecol + yshift


#switch columns (in the case of coordinates are fliped)
swtich_columns <- function(x) {
  tmp <- x$imagecol
  x$imagecol <- x$imagerow
  x$imagerow <- tmp
  return(x)
}
v <- swtich_columns(v)

#save files
saveRDS(v, "")
saveRDS(hires,"")
