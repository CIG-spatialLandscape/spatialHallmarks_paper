OC.st <- readRDS("Desktop/OV_5A.st.rds")
OC.st <- readRDS("Desktop/IJC/datasets/IGTP/4A/RDS/OV_4A.st.rds")
OC.st <- readRDS("Desktop/OV_5B.st.rds")

SaveH5Seurat(st, filename = "Desktop/OV_5A.h5Seurat", overwrite = T)
Convert("Desktop/OV_A.h5Seurat", dest = "h5ad", overwrite = T)
st@images$OV5A@coordinates
st@assays$SCT@counts

write10xCounts("Desktop/matrix.h5", st@assays$SCT@counts)


library(data.table)
data_to_write_out <- as.data.frame(as.matrix(OC.st@assays$Spatial@counts))
fwrite(x = data_to_write_out, row.names = T,file = "matrix_4A_1.csv")
data_to_write_out <- as.data.frame(as.matrix(OC.st@images$OV4A@coordinates))
fwrite(x = data_to_write_out, row.names = T,file = "matrix_coordinates_4A.csv")
m <- read.csv("Desktop/IJC/cell2location/Untitled Folder/data.csv", row.names = 1)
m <- read.csv("Desktop/IJC/cell2location/Untitled Folder/data_4A.csv", row.names = 1)
m <- read.csv("Desktop/IJC/cell2location/Untitled Folder/data_5B.csv", row.names = 1)
m <- read.csv("Desktop/IJC/cell2location/Untitled Folder/data_4A_olbrecht.csv", row.names = 1)
m <- read.csv("Desktop/IJC/cell2location/Untitled Folder/data_public.csv", row.names = 1)
m <- read.csv("Desktop/IJC/cell2location/Untitled Folder/data_public.csv", row.names = 1)
m <- read.csv("Desktop/IJC/cell2location/Untitled Folder/data_4A_integrated_pcg.csv", row.names = 1)
m <- read.csv("Desktop/IJC/cell2location/Untitled Folder/data_4A_integrated_pcg_2.csv", row.names = 1)
OC.st <- AddMetaData(OC.st, m)
colnames(OC.st@meta.data)[11:31]

m <- m/rowSums(m)
m[m < 0.05] <- 0

col.ls <- DiscretePalette(n=21, palette = "polychrome")[3:23]
DimPlot(samples.combined, reduction="umap", cols=col.ls, label=T) + labs(title = "Cell type annotation (Subtype)", caption = "dataset: Olalekan2021")
DimPlot(samples.combined, reduction="umap", cols=col.ls, label=T) + labs(title = "Cell type annotation (Subtype)", caption = "dataset: integration")
DimPlot(samples.combined, reduction="umap", group.by = "MajorCellTypes", cols=col.ls, label=T, repel = T) + labs(title = "Cell type annotation (Majortype)", caption = "dataset: Olalekan2021")
p <- SpatialFeaturePlot(OC.st, colnames(OC.st@meta.data)[11:31], ncol = 3)
pdf(file = "Desktop/IJC/datasets/olbrecht2021/figures/cell2loc_olbrecht_4A_integrated_pcg_2_relative.pdf", width = 15, height = 40)
plot(p)
dev.off()



a <- m %>% group_by(clus) %>% summarise_at(vars(B.cell:pDC), mean)
a <- as.matrix(a[,-1])
rownames(a) <- seq(1,6,1)
heatmap(a)
mdata
order(mdata$annot)


########### Spot light
colnames(OC.st@images$OV4A@coordinates)
plotCorrelationMatrix(as.matrix(m))
plotInteractions(as.matrix(m), "heatmap")
plotInteractions(as.matrix(m), "network")

ct <- colnames(m)

paletteMartin <- c(
  "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", 
  "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", 
  "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")

pal <- colorRampPalette(paletteMartin)(length(ct))
n <- c("blue", "darkorange", "yellow", "purple", "gray", "darkgreen", "lightblue", "darkred")
names(n) <- colnames(m)
names(pal) <- ct
names(col.ls) <- ct

d <- data.frame(row.names  =colnames(OC.st), imagecol=OC.st@images$OV4A@coordinates$imagecol,imagerow= OC.st@images$OV4A@coordinates$imagerow)

df <- merge(m, d, by = 0, all = TRUE)

ggplot() +
  coord_fixed() + scatterpie::geom_scatterpie(
  data = df,
  aes(
    x = imagecol,
    y = -imagerow
  ),
  cols = colnames(m),
  color = NA,
  alpha = 1,
  pie_scale = 0.4)+ 
  theme_void() + 
  theme(legend.key.size = unit(0.5, "lines")) + scale_fill_manual(
    values = col.ls,
    breaks = names(col.ls))


m <- m/rowSums(m)
m <- m[,c(1,2,3,4,9,12,13,14,16,17,18)]
m$Tcells <- m$CD4..Naive.T.cells+m$CD8..T.cells+m$CD8..T.cells...NK+m$T.helper
m <- m[,c(1,5,6,7,8,10,11,12)]
n <- c("#5A5156", "#1CFFCE", "#DEA0FD", "#AA0DFE", "#F8A19F", "#1C8356", "#85660D", "#7EBEFF")
names(n) <- colnames(m)

####
m$Epithelial <- m$Epithelial+m$MKI67..ESCs
m$MKI67..ESCs <- NULL

m$Tcells <- m$CD4..Naive.T.cells+m$CD8..T.cells+m$CD8..T.cells...NK+m$T.helper
m$CD8..T.cells <- NULL
m$CD4..Naive.T.cells <- NULL
m$CD8..T.cells...NK <- NULL
m$T.helper <- NULL

a <- m %>% group_by(clus) %>% summarise_at(vars(B.cell:Tcells), mean)
a <- a %>% group_by(City) %>% mutate(per = Cnt/sum(Cnt))
ggplot(data = a, aes(x = "", y = per, fill = Type)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = Cnt), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  facet_grid(facets=. ~ City)  +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()) + theme(legend.position='bottom') + guides(fill=guide_legend(nrow=2,byrow=TRUE))


###############3
samples.combined2 <- RunUMAP(samples.combined2, return.model = T, dims = 1:30, reduction = "pca")
anchors <- FindTransferAnchors(reference = samples.combined2, query = samples.combined,
                                        dims = 1:30, reference.reduction = "pca", features = rownames(samples.combined))
samples.combined <- MapQuery(anchorset = anchors, reference = samples.combined2, query = samples.combined,
                           refdata = list(celltype = "SubTypes"), reference.reduction = "pca", reduction.model = "umap")
DimPlot(samples.combined2, reduction = "ref.umap", label = T, group.by = "sub.cluster") + 
  DimPlot(samples.combined, reduction = "umap", label = T)

DimPlot(samples.combined2, reduction = "umap", label = T, group.by = "sub.cluster") + 
  DimPlot(samples.combined, reduction = "ref.umap", label = T)

SpatialFeaturePlot(OC.st, features = "Epithelial.Proliferative") + scale_fill_stepsn(n.breaks=6, colours = viridis::viridis(6))
SpatialFeaturePlot(OC.st, features = "Epithelial") + scale_fill_stepsn(n.breaks=6, colours = viridis::viridis(6))
SpatialFeaturePlot(OC.st, features = "TAMs") + scale_fill_stepsn(n.breaks=6, colours = viridis::viridis(6))
