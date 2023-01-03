#####
library(Seurat)
#Get Neighbours


STobject <- readRDS("Desktop/enhanced/new/BreastA_enhanced.rds")


SpatialDimPlot(STobject, pt.size.factor=0.65)

STobject@images[[1]]@coordinates$imagerow

#real distances
coord <- STobject@images[[1]]@coordinates[, c("row", "col")]


for (subspot in rownames(coord)) {
  realrow <-  coord[paste0(substr(subspot, 0, nchar(subspot)-2), ".5"), "row"]
  realrow <- (realrow+1)*100
  realcol <- trunc(coord[paste0(substr(subspot, 0, nchar(subspot)-2), ".3"), "col"])
  realcol <- (realcol+1)*100
  n <- substr(subspot, nchar(subspot), nchar(subspot))
  
  factor <- 55/3
  if (n == 1) {
    realrow <- realrow + factor
    realcol <- realcol + factor
  }
  else if (n == 2) {
    realrow <- realrow + factor
    realcol <- realcol - factor
  }
  else if (n == 3) {
    realrow <- realrow - factor
    realcol <- realcol + factor
  }
  else if (n == 4) {
    realrow <- realrow - factor
    realcol <- realcol - factor
  }
  else if (n == 5) {
    realrow <- realrow
    realcol <- realcol + 2 * factor
  }
  else if (n == 6) {
    realrow <- realrow
    realcol <- realcol - 2 * factor
  }
  coord[subspot, c("realrow", "realcol")] <- c(realrow, realcol)
  
  
}

STobject@images[[1]]@coordinates[, c("realrow", "realcol")] <- coord[, c("realrow", "realcol")]

spot <- 900
STobject$test <- NA
STobject$test[rownames(distances)[spot]] <- "1"
STobject$test[names(distances[spot, distances[spot,] < 200 & distances[spot,] > 0])] <- "2"
SpatialDimPlot(STobject, group.by="test", pt.size.factor=0.65, cols = c("yellow", "blue"))



distances <- as.matrix(dist(coord[, c("realrow", "realcol")], diag = T, upper = T))
names(distances[100, distances[100,] < 120 & distances[100,] > 0])




#####   pv.clus
library(pvclust)



h.clust <- pvclust(t(scale(t(H_NMF_filtered.mt))), method.hclust = "average", parallel = T)

plot(h.clust)
pvrect(h.clust, alpha = 0.95)
print(h.clust)
msplot(h.clust)
pvpick(h.clust, alpha = 0.85)


f.clust <- pvclust(t(H_NMF_filtered.mt), method.hclust = "average", parallel = T)
plot(f.clust)
pvrect(f.clust, alpha = 0.95)
tail(f.clust$edges)

library(cluster)
clu <- diana(t(H_cor_filtered.mt))

clu <- diana(H_cor_filtered.mt)
plot(clu)
pltree(clu)

library(dendextend)
result <- h.clust
result %>%
  as.dendrogram() %>%
  hang.dendrogram() %>%
  plot(main = "Cluster dendrogram with AU/BP values (%)")
result %>% text()
result %>% pvrect(alpha = 0.95)
dend <- as.dendrogram(result)
dend %>%
  pvclust_show_signif_gradient(result, signif_col_fun = colorRampPalette(c("blue", "darkblue", "green"))) %>%
  plot()
