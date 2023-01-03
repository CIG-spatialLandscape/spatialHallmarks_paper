
 hallmark_names <- c("Sustaining Proliferative Signaling", 
                    "Evading Growth Suppressors",
                    "Avoiding Immune Destruction",
                    "Enabling Replicative Immortality",
                    "Tumour-Promoting Inflammation",
                    "Activating Invasion and Metastasis",
                    "Inducing Angiogenesis",
                    "Genome Instability and Mutation",
                    "Resisting Cell Death",
                    "Deregulating Cellular Energetics",
                    "Senescent cells",
                    "Nonmutational Epigenetic reprogramming",
                    "Unlocking Phenotypic Plasticity")



###########################################
ST <- readRDS("Desktop/enhanced/Acinar_enhanced_St_With_Halls.rds")


ST$TME_specific <- "a"

ST$TME_specific[ST$NMF_clusters == "1"] <- "Muscle"
ST$TME_specific[ST$NMF_clusters == "2"] <- "Cancer"
ST$TME_specific[ST$NMF_clusters == "3"] <- "Cancer"
ST$TME_specific[ST$NMF_clusters == "4"] <- "Immune"
ST$TME_specific[ST$NMF_clusters == "5"] <- "Infiltrated"

palette <- c(Tumor="darkred", Immune="darkgreen", ECM="darkorange1",  Mix="darkorchid3", Muscle="red", Infiltrated="yellow")

p <- SpatialDimPlot(ST, group.by="TME_specific", pt.size.factor = 0.65, cols = palette) + theme(legend.position = "none")
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/NMF_clusters/Acinar.pdf")
plot(p)
dev.off()

palette <- RColorBrewer::brewer.pal(5, name = "Paired")

p <- SpatialDimPlot(ST, group.by = "NMF_clusters", cols = palette, pt.size.factor = 0.65) + 
  theme(legend.key = element_blank(), legend.direction = "horizontal", legend.position = "top", legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 4) ))
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/NMF_clusters2/Acinar.pdf")
plot(p)
dev.off()





library(reshape2)
H_modules <- ST@meta.data[,c("NMF_clusters", paste0("H", 1:13))]
H_modules[2:14] <- scale(H_modules[2:14])
mean(H_modules$H12)
sd(H_modules$H12)

H_modules_resh <- melt(H_modules,id.vars='NMF_clusters', measure.vars=paste0("H", 1:13))
colnames(H_modules_resh)[1:2] <- c("cluster", "Hallmark")

p <- ggboxplot(H_modules_resh, x = "cluster", y = "value", ylab = "Hallmark scores", merge = "flip",
               ggtheme = theme_classic(), fill = "Hallmark", palette = "ucscgb",
               legend = "top",
               title = "",
               font.legend = c(15, "plain", "black"),
               font.tickslab = c(15,"black"),
               font.main = c(15, "plain", "black"),
               font.x = c(15, "plain", "black"),
               font.y = c(15, "plain", "black"),
) 

pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/boxplot_NMF/Acinar.pdf")
plot(p)
dev.off()


H_modules <- ST@meta.data[,c("TME_specific", paste0("H", 1:13))]
H_modules[2:14] <- scale(H_modules[2:14])
mean(H_modules$H12)
sd(H_modules$H12)

H_modules_resh <- melt(H_modules,id.vars='TME_specific', measure.vars=paste0("H", 1:13))
colnames(H_modules_resh)[1:2] <- c("cluster", "Hallmark")

p <- ggboxplot(H_modules_resh, x = "cluster", y = "value", ylab = "Hallmark scores", merge = "flip",
               ggtheme = theme_classic(), fill = "Hallmark", palette = "ucscgb",
               legend = "top",
               title = "",
               font.legend = c(15, "plain", "black"),
               font.tickslab = c(15,"black"),
               font.main = c(15, "plain", "black"),
               font.x = c(15, "plain", "black"),
               font.y = c(15, "plain", "black"),
) 
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/boxplot_label/Acinar.pdf")
plot(p)
dev.off()


######################################3

ST <- readRDS("Desktop/enhanced/OV4A_enhanced.rds")


ST$TME_specific <- "a"

ST$TME_specific[ST$NMF_clusters == "1"] <- "Cancer"
ST$TME_specific[ST$NMF_clusters == "2"] <- "Immune"
ST$TME_specific[ST$NMF_clusters == "3"] <- "ECM"
ST$TME_specific[ST$NMF_clusters == "4"] <- "ECM"
ST$TME_specific[ST$NMF_clusters == "5"] <- "ECM"

palette <- c(Cancer="darkred", Immune="darkgreen", ECM="darkorange1",  Mix="darkorchid3", Muscle="red", Infiltrated="yellow")

p <- SpatialDimPlot(ST, group.by="TME_specific", pt.size.factor = 0.65, cols = palette) + theme(legend.position = "none")
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/NMF_clusters/OV4A.pdf")
plot(p)
dev.off()

palette <- RColorBrewer::brewer.pal(12, name = "Paired")[c(6, 8, 11, 1, 4)]

p <- SpatialDimPlot(ST, group.by = "NMF_clusters", cols = palette, pt.size.factor = 0.65)  + 
  scale_fill_manual(labels = c("Cancer", "Immune", "ECM", "ECM (Endothelial)", "ECM"), values = palette) + 
  theme(legend.key = element_blank(), legend.direction = "horizontal", legend.position = "top", legend.title = element_blank(), legend.text=element_text(size=21)) + 
  guides(fill = guide_legend(override.aes = list(size = 4)))
p
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/NMF_clusters2/OV4A.pdf")
plot(p)
dev.off()





library(reshape2)
H_modules <- ST@meta.data[,c("NMF_clusters", paste0("H", 1:13))]
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
               font.legend = c(15, "plain", "black"),
               font.tickslab = c(15,"black"),
               font.main = c(15, "plain", "black"),
               font.x = c(15, "plain", "black"),
               font.y = c(15, "plain", "black"),
)  + scale_fill_manual(labels = c("Cancer", "Immune", "ECM", "ECM (Endothelial)", "ECM"), values = palette) + 

scale_x_discrete(labels=hallmark_names[c(2,4,8,11,12,3,5,9,1,6,7,10,13)]) + theme(axis.text.x =  element_text(angle = -45, hjust = 0), legend.title = element_blank(),
                                                plot.margin = margin(0.5,4.5,0.5,0.5, "cm"))
p2
ST@meta.data[,c(paste0("H", 1:13))] <- scale(ST@meta.data[,c(paste0("H", 1:13))])

SpatialFeaturePlot(ST, features = c("H1"), pt.size.factor = 0.65) + scale_fill_gradientn("", colours = SpatialColors(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 20)) + ggtitle(hallmark_names[12])

for (hallmark in 1:13){
  p <- SpatialFeaturePlot(ST, features = paste0("H", hallmark), pt.size.factor = 0.65) + scale_fill_gradientn("", colours = viridis::inferno(100)) + 
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 20)) + ggtitle(hallmark_names[hallmark])
  pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/figures/H", hallmark, ".pdf"))
  plot(p)
  dev.off()
}


p2
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/boxplot_NMF/OV4A.pdf", width = 15)
plot(p2)
dev.off()

SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))





plot <- p + p2
fig <- ggarrange(p,p2, widths = c(1,1))
annotate_figure(fig, top=text_grob("OV4A", face = "bold"))

H_modules <- ST@meta.data[,c("TME_specific", paste0("H", 1:13))]
H_modules[2:14] <- scale(H_modules[2:14])
mean(H_modules$H12)
sd(H_modules$H12)

H_modules_resh <- melt(H_modules,id.vars='TME_specific', measure.vars=paste0("H", 1:13))
colnames(H_modules_resh)[1:2] <- c("cluster", "Hallmark")

p2 <- ggboxplot(H_modules_resh, x = "cluster", y = "value", ylab = "Hallmark scores", merge = "flip",
               ggtheme = theme_classic(), fill = "Hallmark", palette = "ucscgb",
               legend = "none",
               title = "",
               font.legend = c(15, "plain", "black"),
               font.tickslab = c(15,"black"),
               font.main = c(15, "plain", "black"),
               font.x = c(15, "plain", "black"),
               font.y = c(15, "plain", "black"),
) 
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/boxplot_label/OV4A.pdf")
plot(p)
dev.off()




####################################################################

ST <- readRDS("Desktop/enhanced/NASH_enhanced_St_With_Halls.rds")


ST$TME_specific <- "a"

ST$TME_specific[ST$NMF_clusters == "1"] <- "Tumor"
ST$TME_specific[ST$NMF_clusters == "2"] <- "Tumor"
ST$TME_specific[ST$NMF_clusters == "3"] <- "Mix"
ST$TME_specific[ST$NMF_clusters == "4"] <- "Infiltrated"
ST$TME_specific[ST$NMF_clusters == "5"] <- "Infiltrated"

palette <- c(Tumor="darkred", Immune="darkgreen", ECM="darkorange1",  Mix="darkorchid3", Muscle="red", Infiltrated="yellow")

p <- SpatialDimPlot(ST, group.by="TME_specific", pt.size.factor = 0.65, cols = palette) + theme(legend.position = "none")
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/NMF_clusters/NASH.pdf")
plot(p)
dev.off()

palette <- RColorBrewer::brewer.pal(5, name = "Paired")

p <- SpatialDimPlot(ST, group.by = "NMF_clusters", cols = palette, pt.size.factor = 0.65) + 
  theme(legend.key = element_blank(), legend.direction = "horizontal", legend.position = "top", legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 4) ))
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/NMF_clusters2/NASH.pdf")
plot(p)
dev.off()





library(reshape2)
H_modules <- ST@meta.data[,c("NMF_clusters", paste0("H", 1:13))]
H_modules[2:14] <- scale(H_modules[2:14])
mean(H_modules$H12)
sd(H_modules$H12)

H_modules_resh <- melt(H_modules,id.vars='NMF_clusters', measure.vars=paste0("H", 1:13))
colnames(H_modules_resh)[1:2] <- c("cluster", "Hallmark")

p <- ggboxplot(H_modules_resh, x = "cluster", y = "value", ylab = "Hallmark scores", merge = "flip",
               ggtheme = theme_classic(), fill = "Hallmark", palette = "ucscgb",
               legend = "top",
               title = "",
               font.legend = c(15, "plain", "black"),
               font.tickslab = c(15,"black"),
               font.main = c(15, "plain", "black"),
               font.x = c(15, "plain", "black"),
               font.y = c(15, "plain", "black"),
) 

pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/boxplot_NMF/NASH.pdf")
plot(p)
dev.off()


H_modules <- ST@meta.data[,c("TME_specific", paste0("H", 1:13))]
H_modules[2:14] <- scale(H_modules[2:14])
mean(H_modules$H12)
sd(H_modules$H12)

H_modules_resh <- melt(H_modules,id.vars='TME_specific', measure.vars=paste0("H", 1:13))
colnames(H_modules_resh)[1:2] <- c("cluster", "Hallmark")

p <- ggboxplot(H_modules_resh, x = "cluster", y = "value", ylab = "Hallmark scores", merge = "flip",
               ggtheme = theme_classic(), fill = "Hallmark", palette = "ucscgb",
               legend = "top",
               title = "",
               font.legend = c(15, "plain", "black"),
               font.tickslab = c(15,"black"),
               font.main = c(15, "plain", "black"),
               font.x = c(15, "plain", "black"),
               font.y = c(15, "plain", "black"),
) 
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/boxplot_label/NASH.pdf")
plot(p)
dev.off()


#############################

ST <- readRDS("Desktop/enhanced/HCV1_enhanced_St_With_Halls.rds")


ST$TME_specific <- "a"

ST$TME_specific[ST$NMF_clusters == "1"] <- "Tumor"
ST$TME_specific[ST$NMF_clusters == "2"] <- "Tumor"
ST$TME_specific[ST$NMF_clusters == "3"] <- "Tumor"
ST$TME_specific[ST$NMF_clusters == "4"] <- "Infiltrated"
ST$TME_specific[ST$NMF_clusters == "5"] <- "Infiltrated"

palette <- c(Tumor="darkred", Immune="darkgreen", ECM="darkorange1",  Mix="darkorchid3", Muscle="red", Infiltrated="yellow")

p <- SpatialDimPlot(ST, group.by="TME_specific", pt.size.factor = 0.65, cols = palette) + theme(legend.position = "none")
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/NMF_clusters/HCV1.pdf")
plot(p)
dev.off()

palette <- RColorBrewer::brewer.pal(5, name = "Paired")

p <- SpatialDimPlot(ST, group.by = "NMF_clusters", cols = palette, pt.size.factor = 0.65) + 
  theme(legend.key = element_blank(), legend.direction = "horizontal", legend.position = "top", legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 4) ))
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/NMF_clusters2/HCV1.pdf")
plot(p)
dev.off()





library(reshape2)
H_modules <- ST@meta.data[,c("NMF_clusters", paste0("H", 1:13))]
H_modules[2:14] <- scale(H_modules[2:14])
mean(H_modules$H12)
sd(H_modules$H12)

H_modules_resh <- melt(H_modules,id.vars='NMF_clusters', measure.vars=paste0("H", 1:13))
colnames(H_modules_resh)[1:2] <- c("cluster", "Hallmark")

p <- ggboxplot(H_modules_resh, x = "cluster", y = "value", ylab = "Hallmark scores", merge = "flip",
               ggtheme = theme_classic(), fill = "Hallmark", palette = "ucscgb",
               legend = "top",
               title = "",
               font.legend = c(15, "plain", "black"),
               font.tickslab = c(15,"black"),
               font.main = c(15, "plain", "black"),
               font.x = c(15, "plain", "black"),
               font.y = c(15, "plain", "black"),
) 

pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/boxplot_NMF/HCV1.pdf")
plot(p)
dev.off()


H_modules <- ST@meta.data[,c("TME_specific", paste0("H", 1:13))]
H_modules[2:14] <- scale(H_modules[2:14])
mean(H_modules$H12)
sd(H_modules$H12)

H_modules_resh <- melt(H_modules,id.vars='TME_specific', measure.vars=paste0("H", 1:13))
colnames(H_modules_resh)[1:2] <- c("cluster", "Hallmark")

p <- ggboxplot(H_modules_resh, x = "cluster", y = "value", ylab = "Hallmark scores", merge = "flip",
               ggtheme = theme_classic(), fill = "Hallmark", palette = "ucscgb",
               legend = "top",
               title = "",
               font.legend = c(15, "plain", "black"),
               font.tickslab = c(15,"black"),
               font.main = c(15, "plain", "black"),
               font.x = c(15, "plain", "black"),
               font.y = c(15, "plain", "black"),
) 
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/boxplot_label/HCV1.pdf")
plot(p)
dev.off()



######################3

ST <- readRDS("Desktop/enhanced/HCV2_enhanced_St_With_Halls.rds")


ST$TME_specific <- "a"

ST$TME_specific[ST$NMF_clusters == "1"] <- "Tumor"
ST$TME_specific[ST$NMF_clusters == "2"] <- "ECM"
ST$TME_specific[ST$NMF_clusters == "3"] <- "Infiltrated"
ST$TME_specific[ST$NMF_clusters == "4"] <- "Infiltrated"
ST$TME_specific[ST$NMF_clusters == "5"] <- "Infiltrated"

palette <- c(Tumor="darkred", Immune="darkgreen", ECM="darkorange1",  Mix="darkorchid3", Muscle="red", Infiltrated="yellow")

p <- SpatialDimPlot(ST, group.by="TME_specific", pt.size.factor = 0.65, cols = palette) + theme(legend.position = "none")
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/NMF_clusters/HCV2.pdf")
plot(p)
dev.off()

palette <- RColorBrewer::brewer.pal(5, name = "Paired")

p <- SpatialDimPlot(ST, group.by = "NMF_clusters", cols = palette, pt.size.factor = 0.65) + 
  theme(legend.key = element_blank(), legend.direction = "horizontal", legend.position = "top", legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 4) ))
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/NMF_clusters2/HCV2.pdf")
plot(p)
dev.off()





library(reshape2)
H_modules <- ST@meta.data[,c("NMF_clusters", paste0("H", 1:13))]
H_modules[2:14] <- scale(H_modules[2:14])
mean(H_modules$H12)
sd(H_modules$H12)

H_modules_resh <- melt(H_modules,id.vars='NMF_clusters', measure.vars=paste0("H", 1:13))
colnames(H_modules_resh)[1:2] <- c("cluster", "Hallmark")

p <- ggboxplot(H_modules_resh, x = "cluster", y = "value", ylab = "Hallmark scores", merge = "flip",
               ggtheme = theme_classic(), fill = "Hallmark", palette = "ucscgb",
               legend = "top",
               title = "",
               font.legend = c(15, "plain", "black"),
               font.tickslab = c(15,"black"),
               font.main = c(15, "plain", "black"),
               font.x = c(15, "plain", "black"),
               font.y = c(15, "plain", "black"),
) 

pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/boxplot_NMF/HCV2.pdf")
plot(p)
dev.off()


H_modules <- ST@meta.data[,c("TME_specific", paste0("H", 1:13))]
H_modules[2:14] <- scale(H_modules[2:14])
mean(H_modules$H12)
sd(H_modules$H12)

H_modules_resh <- melt(H_modules,id.vars='TME_specific', measure.vars=paste0("H", 1:13))
colnames(H_modules_resh)[1:2] <- c("cluster", "Hallmark")

p <- ggboxplot(H_modules_resh, x = "cluster", y = "value", ylab = "Hallmark scores", merge = "flip",
               ggtheme = theme_classic(), fill = "Hallmark", palette = "ucscgb",
               legend = "top",
               title = "",
               font.legend = c(15, "plain", "black"),
               font.tickslab = c(15,"black"),
               font.main = c(15, "plain", "black"),
               font.x = c(15, "plain", "black"),
               font.y = c(15, "plain", "black"),
) 
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/boxplot_label/HCV2.pdf")
plot(p)
dev.off()


##################


ST <- readRDS("Desktop/enhanced/HBV_enhanced_St_With_Halls.rds")


ST$TME_specific <- "a"

ST$TME_specific[ST$NMF_clusters == "1"] <- "Tumor"
ST$TME_specific[ST$NMF_clusters == "2"] <- "Immune"
ST$TME_specific[ST$NMF_clusters == "3"] <- "Mix"
ST$TME_specific[ST$NMF_clusters == "4"] <- "Infiltrated"
ST$TME_specific[ST$NMF_clusters == "5"] <- "Infiltrated"

palette <- c(Tumor="darkred", Immune="darkgreen", ECM="darkorange1",  Mix="darkorchid3", Muscle="red", Infiltrated="yellow")

p <- SpatialDimPlot(ST, group.by="TME_specific", pt.size.factor = 0.65, cols = palette) + theme(legend.position = "none")
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/NMF_clusters/HBV.pdf")
plot(p)
dev.off()

palette <- RColorBrewer::brewer.pal(5, name = "Paired")

p <- SpatialDimPlot(ST, group.by = "NMF_clusters", cols = palette, pt.size.factor = 0.65) + 
  theme(legend.key = element_blank(), legend.direction = "horizontal", legend.position = "top", legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 4) ))
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/NMF_clusters2/HBV.pdf")
plot(p)
dev.off()





library(reshape2)
H_modules <- ST@meta.data[,c("NMF_clusters", paste0("H", 1:13))]
H_modules[2:14] <- scale(H_modules[2:14])
mean(H_modules$H12)
sd(H_modules$H12)

H_modules_resh <- melt(H_modules,id.vars='NMF_clusters', measure.vars=paste0("H", 1:13))
colnames(H_modules_resh)[1:2] <- c("cluster", "Hallmark")

p <- ggboxplot(H_modules_resh, x = "cluster", y = "value", ylab = "Hallmark scores", merge = "flip",
               ggtheme = theme_classic(), fill = "Hallmark", palette = "ucscgb",
               legend = "top",
               title = "",
               font.legend = c(15, "plain", "black"),
               font.tickslab = c(15,"black"),
               font.main = c(15, "plain", "black"),
               font.x = c(15, "plain", "black"),
               font.y = c(15, "plain", "black"),
) 

pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/boxplot_NMF/HBV.pdf")
plot(p)
dev.off()


H_modules <- ST@meta.data[,c("TME_specific", paste0("H", 1:13))]
H_modules[2:14] <- scale(H_modules[2:14])
mean(H_modules$H12)
sd(H_modules$H12)

H_modules_resh <- melt(H_modules,id.vars='TME_specific', measure.vars=paste0("H", 1:13))
colnames(H_modules_resh)[1:2] <- c("cluster", "Hallmark")

p <- ggboxplot(H_modules_resh, x = "cluster", y = "value", ylab = "Hallmark scores", merge = "flip",
               ggtheme = theme_classic(), fill = "Hallmark", palette = "ucscgb",
               legend = "top",
               title = "",
               font.legend = c(15, "plain", "black"),
               font.tickslab = c(15,"black"),
               font.main = c(15, "plain", "black"),
               font.x = c(15, "plain", "black"),
               font.y = c(15, "plain", "black"),
) 
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/boxplot_label/HBV.pdf")
plot(p)
dev.off()

