
files <- list.files(path="Desktop/enhanced", pattern="*.rds", full.names=TRUE, recursive=FALSE)
file <- files[1]
STobject <- readRDS(file)
name <- strsplit(strsplit(file, "/")[[1]][3], "_")[[1]][1] 
name
gc()
labels <- c("Muscle", "Cancer", "Cancer", "Immune", "Infilitrated")
palette <- RColorBrewer::brewer.pal(12, name = "Paired")[c(6, 8, 11, 1, 4)]

p <- SpatialDimPlot(STobject, group.by = "NMF_clusters", cols = palette, pt.size.factor = 0.65)  + 
  scale_fill_manual(labels = labels, values = palette) + 
  theme(legend.key = element_blank(), legend.direction = "horizontal", legend.position = "top", legend.title = element_blank(), legend.text=element_text(size=30)) + 
  guides(fill = guide_legend(override.aes = list(size = 4)))
p
pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/figures/NMF_clusters2/", name, ".pdf"), width = 10, height = 7)
plot(p)
dev.off()

library(reshape2)
H_modules <- STobject@meta.data[,c("NMF_clusters", paste0("H", 1:13))]
H_modules[2:14] <- scale(H_modules[2:14])
mean(H_modules$H12)
sd(H_modules$H12)

H_modules_resh <- melt(H_modules,id.vars='NMF_clusters', measure.vars=paste0("H", 1:13))
colnames(H_modules_resh)[1:2] <- c("cluster", "Hallmark")
H_modules_resh$cluster <- as.factor(H_modules_resh$cluster)
H_modules_resh$Hallmark <- factor(H_modules_resh$Hallmark, 
                                  levels=c("H2", "H4", "H8", "H11", "H12", "H3", "H5", "H9", "H1", "H6", "H7", "H10", "H13"))

hallmark_names[c(2,4,8,11,12,3,5,9,1,6,7,10,13)]

p2 <- ggboxplot(H_modules_resh, x = "Hallmark", y = "value", ylab = "Scaled Hallmark scores", merge = "flip",
                ggtheme = theme_classic(), fill = "cluster",
                legend = "top",
                title = "", xlab = "", outlier.shape = 20,
                font.legend = c(20, "plain", "black"),
                font.tickslab = c(15,"black"),
                font.main = c(15, "plain", "black"),
                font.x = c(35, "plain", "black"),
                font.y = c(15, "plain", "black"),
)  + scale_fill_manual(labels = labels, values = palette) + 
  
  scale_x_discrete(labels=hallmark_names[c(2,4,8,11,12,3,5,9,1,6,7,10,13)]) + theme(axis.text.x =  element_text(angle = -45, hjust = 0), legend.title = element_blank(),
                                                                                    plot.margin = margin(0.5,4.5,0.5,0.5, "cm"))
p2

pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/figures/boxplot_NMF/", name, ".pdf"), width = 10, height = 8)
plot(p2)
dev.off()


STobject@meta.data[,paste0("factor_",1:5)] <- scale(STobject@meta.data[,paste0("factor_",1:5)])
for (i in 1:5) {
  p <- SpatialFeaturePlot(STobject, features = paste0("factor_",i), pt.size.factor = 0.65) + scale_fill_gradientn("", colours = viridis::inferno(100)) + 
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 40)) + ggtitle(paste0("factor_",i))
  pdf(paste0("Desktop/IJC/datasets/Public/",name, "/figures/factors5/factor_", i, ".pdf"), width = 10, height = 7)
  plot(p)
  dev.off()
}

STobject@meta.data[,paste0("H",1:13)] <- scale(STobject@meta.data[,paste0("H",1:13)])
##
for (i in 1:13) {
  p <- SpatialFeaturePlot(STobject, features = paste0("H",i), pt.size.factor = 0.65) + scale_fill_gradientn("", colours = viridis::inferno(100)) + 
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 30)) + ggtitle(hallmark_names[i])
  pdf(paste0("Desktop/finalplots/H", i, ".pdf"), width = 10, height = 7)
  plot(p)
  dev.off()
}


p
