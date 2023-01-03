OC$f15 <- "0"
factor <- OC@reductions$NMF@cell.embeddings[,15]
x <- mean(factor) + 1.5*sd(factor)
OC$f15[OC@reductions$NMF@cell.embeddings[,15] > x] <- "1"
SpatialDimPlot(OC, group.by = c("f15", "bayes.cluster"))


for (i in 1:15) {
  OC[[paste0("f",i,"_label")]][1] <- "0"
  factor <- OC@reductions$NMF@cell.embeddings[,i]
  x <- mean(factor) + 1.5*sd(factor)
  OC[[paste0("f",i,"_label")]][OC@reductions$NMF@cell.embeddings[,i] > x,] <- "1"
}
