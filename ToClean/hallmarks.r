H <- read.table("Desktop/Hallmarks_genes_PthCommons2019_final_onlyGenes.txt", sep = "\t", header = T)


H_features <- list("1"= H[H$H =="H1",1],"2"= H[H$H =="H2",1], "3"= H[H$H =="H3",1],"4"=H[H$H =="H4",1],"5"= H[H$H =="H5",1],"6"= H[H$H =="H6",1],"7"= H[H$H =="H7",1],
                   "8"= H[H$H =="H8",1],"9"= H[H$H =="H9",1],"10"= H[H$H =="H10",1],"11"= H[H$H =="H11",1],"12"= H[H$H =="H12",1], "13"= H[H$H =="H13",1])


OV_4A_enhanced.St <- AddModuleScore(
  OV_4A_enhanced.St,
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

path_figures <- "Desktop/IJC/datasets/IGTP/4A/figures/"
colnames(OV_4A_enhanced.St@meta.data)[19:32]
for (h in colnames(OV_4A_enhanced.St@meta.data)[19:32]){
  file_name = paste0(path_figures, "hallmarks/", h, ".png")
  ggsave(file_name, plot = SpatialFeaturePlot(OV_4A_enhanced.St, features = h, pt.size.factor = 0.6), height = 7, width = 7, bg = "white")
}



colnames(OV_4A_enhanced.St@meta.data)
H_modules <- OV_4A_enhanced.St@meta.data[,18:31]
H_modules[2:14] <- scale(H_modules[2:14])
mean(H_modules$H12)
sd(H_modules$H12)

library(reshape2)
H_modules_resh <- melt(H_modules,id.vars='NMF_clusters', measure.vars=colnames(H_modules[2:14]))
colnames(H_modules_resh)[1:2] <- c("NMF_clusters", "Hallmark")

ggboxplot(H_modules_resh, x = "cluster", y = "value", ylab = "Hallmark scores", merge = "flip",
          ggtheme = theme_classic(), fill = "Hallmark", palette = "ucscgb",
          ylim = c(-4,7), legend = "top",
          title = "",
          font.legend = c(15, "plain", "black"),
          font.tickslab = c(15,"black"),
          font.main = c(15, "plain", "black"),
          font.x = c(15, "plain", "black"),
          font.y = c(15, "plain", "black"),
) 

ggplot(H_modules_resh) + geom_boxplot(aes(x=Hallmark, y = value, fill = Hallmark)) + facet_wrap(~NMF_clusters)
