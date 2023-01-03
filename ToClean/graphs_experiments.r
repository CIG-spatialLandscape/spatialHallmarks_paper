library(Seurat)


STobject <- readRDS("Desktop/IJC/datasets/Dutreneo/RDS/Assays/DU2.rds")
STobject
DefaultAssay(STobject) <- "Deconvolution"
SpatialFeaturePlot(STobject, features = "Epithelium")
SpatialFeaturePlot(STobject, features = "Fibroblast")
SpatialFeaturePlot(STobject, features = "Tcell")



Tcell <- kmeans(t(as.matrix(STobject@assays$Deconvolution["Tcell",])), centers = 2)$cluster
STobject$Tcell <- Tcell
SpatialDimPlot(STobject, group.by = "Tcell")

t(as.matrix(STobject@assays$Deconvolution["Tcell",])) > 0.1 & t(as.matrix(STobject@assays$Deconvolution["Epithelium",])) > 0.3
shared <- t(as.matrix(STobject@assays$Deconvolution["Tcell",])) > 0.1 & t(as.matrix(STobject@assays$Deconvolution["Epithelium",])) > 0.5
shared <- t(as.matrix(STobject@assays$Deconvolution["Fibroblast",])) > 0.2 & t(as.matrix(STobject@assays$Deconvolution["Epithelium",])) > 0.6

STobject$shared <- "1"
STobject$shared[shared] <- "2"
SpatialDimPlot(STobject, group.by = "shared")


library(ggplot2)
df <- t(as.matrix(STobject@assays$Deconvolution@data))
library(reshape2)
df1 <- melt(df)
df1$range <- 1:nrow(df1)
ggplot(df1, aes(x=Var1, y=value, col=Var2)) + geom_line()






library(igraph)
d <- as.matrix(dist(STobject@images$DU2@coordinates[,c("row", "col")], diag = T, upper = T))

d[d > 3] <- 0


g <- graph_from_adjacency_matrix(d, weighted=TRUE, mode = "undirected")
g <- set.vertex.attribute(g, "x", value=STobject@images$DU2@coordinates$imagecol) 
g <- set.vertex.attribute(g, "y", value=-STobject@images$DU2@coordinates$imagerow)

plot(g, vertex.size = 1, vertex.label.cex = 0.001)



# t cells and epithelial connexion
shared <- t(as.matrix(STobject@assays$Deconvolution["Tcell",])) > 0.1 & t(as.matrix(STobject@assays$Deconvolution["Epithelium",])) > 0.5
d[!shared, !shared] <- 0

g <- graph_from_adjacency_matrix(d, weighted=TRUE, mode = "undirected")
g <- set.vertex.attribute(g, "x", value=STobject@images$DU2@coordinates$imagecol) 
g <- set.vertex.attribute(g, "y", value=-STobject@images$DU2@coordinates$imagerow)

plot(g, vertex.size = 1, vertex.label.cex = 0.001)


component_distribution(g)

 c <- components(g, mode = c("weak", "strong"))

is_connected(g, mode = c("weak", "strong"))

count_components(g, mode = c("weak", "strong"))


rows_to_keep <- rowSums(d) != 0
cols_to_keep <- colSums(d) != 0

d2 <- d[shared, shared]
g <- graph_from_adjacency_matrix(d2, weighted=TRUE, mode = "undirected")
g <- set.vertex.attribute(g, "x", value=STobject@images$DU2@coordinates$imagecol[shared]) 
g <- set.vertex.attribute(g, "y", value=-STobject@images$DU2@coordinates$imagerow[shared])

plot(g, vertex.size = 1, vertex.label.cex = 0.001)

which.max(c$csize)
c$membership==10

names(c$membership)[c$membership==10]

STobject$test <- NA
STobject$test[names(c$membership)[c$membership==20]] <- "Component1"
STobject$test[names(c$membership)[c$membership==23]] <- "Component2"
STobject$test[names(c$membership)[c$membership!=20 & c$membership!=23]] <- "Component3"
SpatialDimPlot(STobject, group.by = "test")
Idents(STobject) <- STobject$test
DefaultAssay(STobject) <- "SCT"
m <- FindAllMarkers(STobject, only.pos = T)

SpatialFeaturePlot(STobject, features = "SH3GL2")


d1 <- d2
d1[d1 > 0] <- 1
hist(rowSums(d1))
degree <- as.data.frame(rowSums(d1))
hist(degree$`rowSums(d1)`)

DefaultAssay(STobject) <- "dorothea"
m <- FindAllMarkers(STobject, only.pos = T, test.use = "roc")
SpatialFeaturePlot(STobject, features =  "MYB")


STobject$test1 <- NA
STobject$test1[names(c$membership)[c$membership==20]] <- "Component"
STobject$test1[rownames(degree)[degree>8 & rownames(degree) %in% names(c$membership)[c$membership==20]]] <- "Component-HD"
SpatialDimPlot(STobject, group.by = "test1")



epithelial <- t(as.matrix(STobject@assays$Deconvolution["Epithelium",])) > 0.4
STobject$epithelial <- "0"
STobject$epithelial[epithelial] <- "1"
SpatialDimPlot(STobject, group.by = "epithelial")


Fibroblast <- t(as.matrix(STobject@assays$Deconvolution["Fibroblast",])) > 0.4
STobject$Fibroblast <- "0"
STobject$Fibroblast[Fibroblast] <- "1"
SpatialDimPlot(STobject, group.by = "Fibroblast")



STobject$experiment <- NA
STobject$experiment[epithelial] <- "Epithelial"
STobject$experiment[Fibroblast] <- "Fibroblast"
STobject$experiment[names(c$membership)[c$membership==20]] <- "Tcell infilitration"
SpatialDimPlot(STobject, group.by = "experiment")

nodes <- unique(c(colnames(STobject)[epithelial], colnames(STobject)[epithelial]))

d <- as.matrix(dist(STobject@images$DU2@coordinates[,c("row", "col")], diag = T, upper = T))

#d <- d[nodes, nodes]
d[colnames(STobject)[epithelial], colnames(STobject)[epithelial]] <- 0
d[colnames(STobject)[Fibroblast], colnames(STobject)[Fibroblast]] <- 0
d[d > 3] <- 0

rows_to_keep <- rowSums(d) != 0
cols_to_keep <- colSums(d) != 0
d <- d[rows_to_keep, cols_to_keep] 

g <- graph_from_adjacency_matrix(d, weighted=TRUE, mode = "undirected")
g <- set.vertex.attribute(g, "x", value=STobject@images$DU2@coordinates$imagecol[rows_to_keep]) 
g <- set.vertex.attribute(g, "y", value=-STobject@images$DU2@coordinates$imagerow[cols_to_keep])

plot(g, vertex.size = 1, vertex.label.cex = 0.001)



DU14 <- readRDS("Desktop/IJC/datasets/Dutreneo/RDS/Assays/DU14.rds")

epithelial <- t(as.matrix(DU14@assays$Deconvolution["Epithelium",])) > 0.4
DU14$epithelial <- "0"
DU14$epithelial[epithelial] <- "1"
SpatialDimPlot(DU14, group.by = "epithelial")


Fibroblast <- t(as.matrix(DU14@assays$Deconvolution["Fibroblast",])) > 0.4
DU14$Fibroblast <- "0"
DU14$Fibroblast[Fibroblast] <- "1"
SpatialDimPlot(DU14, group.by = "Fibroblast")



DU14$experiment <- NA
DU14$experiment[epithelial] <- "Epithelial"
DU14$experiment[Fibroblast] <- "Fibroblast"
DU14$experiment[names(c$membership)[c$membership==20]] <- "Tcell infilitration"
SpatialDimPlot(DU14, group.by = "experiment")

nodes <- unique(c(colnames(DU14)[epithelial], colnames(DU14)[epithelial]))

d <- as.matrix(dist(DU14@images$DU14@coordinates[,c("row", "col")], diag = T, upper = T))

#d <- d[nodes, nodes]
d[colnames(DU14)[epithelial], colnames(DU14)[epithelial]] <- 0
d[colnames(DU14)[Fibroblast], colnames(DU14)[Fibroblast]] <- 0
d[d > 3] <- 0

rows_to_keep <- rowSums(d) != 0
cols_to_keep <- colSums(d) != 0
d <- d[rows_to_keep, cols_to_keep] 

g2 <- graph_from_adjacency_matrix(d, weighted=TRUE, mode = "undirected")
g2 <- set.vertex.attribute(g2, "x", value=DU14@images$DU14@coordinates$imagecol[rows_to_keep]) 
g2 <- set.vertex.attribute(g2, "y", value=-DU14@images$DU14@coordinates$imagerow[cols_to_keep])

plot(g2, vertex.size = 1, vertex.label.cex = 0.001)


components(g2)
d1 <- d
d1[d1 > 0] <- 1
hist(rowSums(d1))
degree <- as.data.frame(rowSums(d1))

group1 <- rownames(degree)[degree == 1]
group2 <- rownames(degree)[degree > 1 & degree < 4]
group3 <- rownames(degree)[degree >= 4 & degree < 6]
group4 <- rownames(degree)[degree >= 6 & degree < 8]
group5 <- rownames(degree)[degree >= 8 & degree < 12]
group6 <- rownames(degree)[degree == 12]

DU14$test <- NA
DU14@meta.data[group1, ]$test <- "1"
DU14@meta.data[group2, ]$test <- "2-3"
DU14@meta.data[group3, ]$test  <- "4-5"
DU14@meta.data[group4, ]$test <- "6-7"
DU14@meta.data[group5, ]$test <- "8-11"
DU14@meta.data[group6, ]$test  <- "12"
SpatialDimPlot(DU14, group.by = "test")

VlnPlot(DU14, features = paste0("H", 1:13), group.by = "test")
SpatialFeaturePlot(DU14, features = "Tcell")
