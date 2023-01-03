OC$f15 <- "0"
x <- quantile(OC@reductions$NMF@cell.embeddings[,15], 0.9)
x <- quantile(OC@reductions$NMF@cell.embeddings[,15], 0.95)
OC$f15[OC@reductions$NMF@cell.embeddings[,15] > x] <- "1"
SpatialDimPlot(OC, group.by = "f15")

OC$f1 <- "0"
x <- quantile(OC@reductions$NMF@cell.embeddings[,1], 0.9)
x <- quantile(OC@reductions$NMF@cell.embeddings[,1], 0.95)
OC$f1[OC@reductions$NMF@cell.embeddings[,1] > x] <- "1"
SpatialDimPlot(OC, group.by = "f1")

OC$f9 <- "0"
x <- quantile(OC@reductions$NMF@cell.embeddings[,5], 0.90)
OC$f9[OC@reductions$NMF@cell.embeddings[,5] > x] <- "1"
SpatialDimPlot(OC, group.by = "f9")

OC$idents <- "0"
OC$idents[OC$f1 == "1"] <- "f1"
OC$idents[OC$f15 == "1"] <- "f15"
OC$idents[OC$f9 == "1"] <- "f9"
SpatialDimPlot(OC, group.by = "idents")

OC$f2 <- "0"
x <- quantile(OC@reductions$NMF@cell.embeddings[,2], 0.9)
x <- quantile(OC@reductions$NMF@cell.embeddings[,2], 0.95)
OC$f2[OC@reductions$NMF@cell.embeddings[,2] > x] <- "1"
SpatialDimPlot(OC, group.by = "f2")

Idents(OC) <- OC$idents
m <- FindMarkers(OC, only.pos = T, ident.1 = "1", logfc.threshold = 0.5)
m <- FindAllMarkers(OC, only.pos = T)
m <- m[order(m$avg_log2FC, decreasing = T),]
df1 <- gost(query = rownames(m), organism = "hsapiens", sources = "GO:BP", evcodes = TRUE, ordered_query = T )
df2 <- gost(query = rownames(m), organism = 'gp__aJd5_EBf2_fFc', evcodes = TRUE)$result
gostplot(df1)
g <- m[m$cluster=="f15",]$gene
df1 <- gost(query = rownames(m), organism = "hsapiens", sources = "GO:BP", evcodes = TRUE)

m <- FindMarkers(OC, only.pos = T, ident.1 = "f15", ident.2 = c("f9"))


DefaultAssay(OC) <- "prediction"

coord <- OC@images$OC@coordinates[,4:5]
dis <- as.matrix(dist(coord, diag = T, upper = T))
pred <- OC@assays$prediction@data
pred <- pred[-nrow(pred),]
diag(dis) <- 0.1
pred[1,1:5]

sum(pred[5,]/dis[5,])

x <- c()
for (i in 1:ncol(pred)){
  x <- c(x, pred[16, i]*sum(pred[5,]/dis[5,]))
}

x1 <- c()
for (i in 1:ncol(pred)){
  x1 <- c(x1, pred[16, i]*sum(pred[9,]/dis[9,]))
}

mat <- matrix(0, nrow = nrow(pred), ncol = nrow(pred))
for (spot in 1:ncol(pred)){
  print(spot)
  for (celltype1 in 1:nrow(mat)) {
    for (celltype2 in celltype1:nrow(mat)) {
      mat[celltype1, celltype2] <- mat[celltype1, celltype2]+pred[celltype1, spot]*sum(pred[celltype2,]/dis[celltype2,])
    }
  }
}

for (spot in 1:ncol(pred)){
  print(spot)
  for (celltype1 in 1:nrow(pred)) {
    for (celltype2 in celltype1:nrow(pred)) {
      l[[celltype1]][[celltype2]] <- c(l[[celltype1]][[celltype2]], pred[celltype1, spot]*sum(pred[celltype2,]/dis[celltype2,]))
    }
  }
}

rownames(mat) <- rownames(pred)
colnames(mat) <- rownames(pred)
mat[lower.tri(mat)] = t(mat)[lower.tri(mat)]
mat[lower.tri(mat)] = 0


l <- list()
for (celltype in rownames(pred)) {
  l[[celltype]] <- list()
  for (celltype2 in rownames(pred)) {
    l[[celltype]][[celltype2]] <- c(0)
  }
}

l <- list(rownames(pred))

m <- matrix(, nrow = 2, ncol = 2)

thresh <- quantile(l$CAFs$CAFs, 0.75)
sum(l$CAFs$CAFs[l$CAFs$CAFs>thresh])

thresh <- quantile(l$Epithelial$Epithelial, 0.99)
sum(l$Epithelial$Epithelial[l$Epithelial$Epithelial>thresh])


thresh <- quantile(l$`CAFs-Immune2`$Epithelial, 0.95)
sum(l$`CAFs-Immune2`$Epithelial[l$`CAFs-Immune2`$Epithelial>thresh])

thresh <- quantile(l$`CAFs-Immune2`$`CAFs-Immune2`, 0.99)
sum(l$`CAFs-Immune2`$`CAFs-Immune2`[l$`CAFs-Immune2`$`CAFs-Immune2`>thresh])


thresh <- quantile(l$`CAFs-Immune2`$Epithelial, 0.99)
sum(l$`CAFs-Immune2`$Epithelial[l$`CAFs-Immune2`$Epithelial>thresh])


d <- l$`CAFs-Immune2`$Epithelial
thresh <- mean(d)+3*sd(d)
sum(d[d>thresh])

embeddings <- scale(OC@reductions$NMF@cell.embeddings)
list_factors <- c()
for (factor in 1:15) {
  
  f <- embeddings[,factor]
  threshold <- mean(f) + 1.5*sd(f)
  top.spots <- colnames(OC)[f > threshold]
  clusters.spots <- colnames(OC)[OC$bayes.cluster %in% c(1,2,5)]
  
  t <- table(top.spots %in% clusters.spots)
  if (t[2]/sum(t) > 0.7 || is.na(t[2]/sum(t))) {
    list_factors <- c(list_factors, factor)
  }

}


OC$f15 <- "0"
factor <- OC@reductions$NMF@cell.embeddings[,15]
x <- mean(factor) + 1.5*sd(factor)
OC$f15[OC@reductions$NMF@cell.embeddings[,15] > x] <- "1"
SpatialDimPlot(OC, group.by = c("f15", "bayes.cluster"))


OC$f1 <- "0"
factor <- OC@reductions$NMF@cell.embeddings[,1]
x <- mean(factor) + 1.5*sd(factor)
OC$f1[OC@reductions$NMF@cell.embeddings[,1] > x] <- "1"
SpatialDimPlot(OC, group.by = c( "f1"))


pred <- t(embeddings[,list_factors]


mat <- matrix(0, nrow = nrow(pred), ncol = nrow(pred))
for (spot in 1:ncol(pred)){
  print(spot)
  for (celltype1 in 1:nrow(mat)) {
    for (celltype2 in celltype1:nrow(mat)) {
      mat[celltype1, celltype2] <- mat[celltype1, celltype2]+pred[celltype1, spot]*sum(pred[celltype2,]/dis[celltype2,])
    }
  }
}



##########

OC <- OC[,colnames(OC)[OC$bayes.cluster %in% c(1,2,5)]]
OC <- SCTransform(OC,return.only.var.genes = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.1, assay = "Spatial")


OC$overlap <- "0"
OC$overlap[OC$f15 == "1" & OC$sub.cluster=="4"] <- "1"

OC$overlap2 <- "0"
OC$overlap2[OC$f1 == "1" & OC$overlap=="1"] <- "1"



###########

OC$overlap3 <- "0"
OC$overlap3[OC$f1 == "1" & OC$sub.cluster=="11"] <- "1"




OC$expriment <- "0"
OC$expriment[OC$f1 == "1"] <- "factor1"
OC$expriment[OC$overlap == "1"] <- "bs_f15"
OC$expriment[OC$overlap2 == "1"] <- "f1_bs_f15"


OC$experiment <- "0"
OC$experiment[OC$overlap == "1"] <- "bs_f15"
OC$experiment[OC$overlap3 == "1"] <- "bs_f1"

stringr::str_replace(colnames(ST), pattern = "_1", "")

OC$nb <- ST
OC$nb_experiment <- OC$nb
OC$nb_experiment[OC$nb_experiment == "nbs_bs_f15" & OC$sub.cluster == "14"] <- "14"
OC$nb_experiment[OC$nb_experiment == "nbs_bs_f15" & OC$sub.cluster == "1"] <- "1"
OC$nb_experiment[OC$nb_experiment == "nbs_bs_f15" & OC$sub.cluster == "3"] <- "3"
OC$nb_experiment[OC$nb_experiment == "nbs_bs_f15" & OC$sub.cluster == "4"] <- "4"
OC$nb_experiment[OC$nb_experiment == "nbs_bs_f15" & OC$sub.cluster == "8"] <- "8"
OC$nb_experiment[OC$nb_experiment == "nbs_bs_f15" & OC$sub.cluster == "11"] <- "11"
OC$nb_experiment[is.na(OC$nb_experiment)] <- "0"
m <- FindMarkers(OC, ident.1 = "bs_f15", ident.2 = "14")
