
OC <- RunPCA(OC, assay = "SCT")
OC <- FindNeighbors(OC, reduction = "pca")
OC <- FindClusters(OC, graph.name = "SCT_snn", resolution = 0.3)
SpatialDimPlot(OC)

library(NICHES)
STobject  <- readRDS("Desktop/IJC/TFG/5merged.rds")



STobject@meta.data$x <- STobject@images[[1]]@coordinates$row
STobject@meta.data$y <- STobject@images[[1]]@coordinates$col

STobject <- SeuratWrappers::RunALRA(STobject)


NICHES_output <- RunNICHES(object = STobject,
                           LR.database = "omnipath",
                           species = "human",
                           assay = "alra",
                           position.x = 'x',
                           position.y = 'y',
                           rad.set = 2, # Geometry dependent
                           min.cells.per.ident = 0,
                           min.cells.per.gene = NULL,
                           meta.data.to.map = c('spatial.cluster'),
                           CellToCell = F,CellToSystem = F,SystemToCell = F,
                           CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = T)

niche <- NICHES_output[['NeighborhoodToCell']]
Idents(niche) <- niche[['ReceivingType']]

# Scale and visualize
niche <- ScaleData(niche)
niche <- FindVariableFeatures(niche,selection.method = "disp")
niche <- RunPCA(niche)
ElbowPlot(niche,ndims = 50)
niche <- RunUMAP(niche,dims = 1:10)
DimPlot(niche,reduction = 'umap',pt.size = 0.5,shuffle = T, label = T) +ggtitle('Cellular Microenvironment')+NoLegend()
mark <- FindAllMarkers(niche,min.pct = 0.25,only.pos = T,test.use = "roc")
GOI_niche <- mark %>% group_by(cluster) %>% top_n(10,myAUC)
DoHeatmap(STobject,features = unique(GOI_niche$gene))+ 
  scale_fill_gradientn(colors = c("grey","white", "blue")) 


# Add Niches output as an assay
niches.data <- GetAssayData(object =  niche[['NeighborhoodToCell']], slot = 'data')
colnames(niches.data) <- niche[['ReceivingCell']]$ReceivingCell
STobject[["NeighborhoodToCell"]] <- CreateAssayObject(data = niches.data )
DefaultAssay(STobject) <- "NeighborhoodToCell"
STobject <- ScaleData(STobject)


SpatialFeaturePlot(STobject,
                   features = c('COL3A1—ITGA1-ITGB1', 'COL6A3—ITGA1-ITGB1'),
                   slot = 'scale.data')
DefaultAssay(STobject) <- "Deconvolution"

SpatialFeaturePlot(STobject, features = c("Fibroblast", 'COL6A3—ITGA1-ITGB1'))
SpatialFeaturePlot(STobject, features = c("Epithelium", 'PSAP—GPR37L1'))
mark[mark$cluster=="11",]$gene[1:5] 
SpatialFeaturePlot(STobject,
                   features = mark[mark$cluster=="2",]$gene[1:5] ,
                   slot = 'scale.data')


STobject$adjacent <- "0"
STobject$adjacent[STobject$Cancer=="Cancer (TME-adjacent)" | STobject$TME=="TME (Cancer-adjacent)"] <- "interface"
SpatialDimPlot(STobject, group.by = "adjacent", pt.size.factor = 0.65)

Idents(STobject) <- STobject$clus
mark2 <- FindAllMarkers(STobject,min.pct = 0.25,only.pos = T,test.use = "roc", assay = "NeighborhoodToCell")


View(mark)
SpatialFeaturePlot(STobject,
                   features = GOI_niche[GOI_niche$cluster!=0,]$gene[6:10] ,
                   slot = 'scale.data', pt.size.factor = 0.65)

VlnPlot(STobject, features = "COL1A2—DDR1", group.by = "clus")



STobject$clus <- "0"
STobject$clus[STobject$Cancer=="Cancer (TME-adjacent)"] <- "Interface"
STobject$clus[STobject$Cancer=="Cancer (TME-distant)"] <- "Cancer (TME-distant)"
STobject$clus[STobject$TME=="TME (Cancer-adjacent)"] <- "Interface"
STobject$clus[STobject$TME=="TME (Cancer-distant)"] <- "TME (Cancer-distant)"
SpatialDimPlot(STobject, group.by = "clus", pt.size.factor = 0.65)
Idents(STobject) <- STobject$clus


DefaultAssay(STobject) <- "SCT"
VlnPlot(STobject, features = c("MSN", "ICAM1"), group.by = "clus")
VlnPlot(STobject, features = "ICAM1", group.by = "clus")

interface_spots <- colnames(STobject)[STobject$clus=="Interface"]

list_correlations <- lapply(paste0("H",1:13), function(h) {
  sort(sapply(rownames(STobject@assays$NeighborhoodToCell), function(x){
    cor(STobject@assays$NeighborhoodToCell@scale.data[x,], STobject@meta.data[,h])
  }), decreasing = T)
})


top <- lapply(list_correlations, function(x){
  head(names(sort(x, decreasing = T)),5)
})



VlnPlot(STobject, features = top[[3]], group.by = "clus")




intersect(top)
Reduce(intersect, top[c(13,11)])


Reduce(intersect, top[c(2,4,8,9,11,12)]) #cancer
Reduce(intersect, top[c(1,3,5,6,7,10,13)]) #TME


common_cancer <- Reduce(intersect, top[c(2,4,8,11,12)]) #cancer
list_correlations[c(2,4,8,11,12)]
common_cor_cancer <- lapply(c(2,4,8,11,12), function(x){
  sort(list_correlations[[x]][common_cancer], decreasing = T)
})
names(common_cor_cancer) <- paste0("H", c(2,4,8,11,12))

Reduce(intersect, top[c(1,5,6,7,10)]) #TME

venn.diagram(top, category.names = paste0("H", 1:13), filename = "Desktop/venn.png" )
