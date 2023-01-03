library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
#summary_comparisions <- data.frame(row.names = 1)
#summary_comparisions$Cancer <- "a"
#summary_comparisions$Cancer_factor <- "a"
#summary_comparisions$Hallmark <- "a"
#summary_comparisions$Diff <- "a"
#summary_comparisions$Sgnif <- "a"

rm(STobject)
files <- list.files(path="Desktop/enhanced", pattern="*.rds", full.names=TRUE, recursive=FALSE)
file <- files[1]
STobject <- readRDS(file)
gc()
factor <- STobject@meta.data[,paste0("factor_", 1:5)]

name <- strsplit(strsplit(file, "/")[[1]][3], "_")[[1]][1] 

name
factor_spots <- lapply(1:5, function(x) {
  colnames(STobject)[factor[,x] > mean(factor[,x]) +1*sd(factor[,x])]
})

cancer_factors <- c(1)
TME_factors <- c(3,4)


cancer <- sapply(cancer_factors, function(x){
  factor_spots[[x]]
})
cancer <- do.call(c, cancer)
cancer <- cancer[!duplicated(cancer)]

TME <- sapply(TME_factors, function(x){
  factor_spots[[x]]
})
TME <- do.call(c, TME)
TME <- TME[!duplicated(TME)]

duplicated <- c(TME, cancer)[duplicated(c(TME, cancer))]
TME <- TME[!TME %in% duplicated]
cancer <- cancer[!cancer %in% duplicated]



coord <- STobject@images[[1]]@coordinates[c(cancer,TME),]
distances <- as.matrix(dist(coord, diag = T, upper = T))
distances <- distances[cancer, TME]

# Cancer cells
min_dist <- apply(distances,1,min)

inner <- names(min_dist[min_dist > quantile(min_dist, 0.75)])
outer <- names(min_dist[min_dist < quantile(min_dist, 0.25)])


STobject$Cancer <- NA
STobject$Cancer[inner] <- "TME-distant"
STobject$Cancer[outer] <- "TME-proximal"

p <- SpatialDimPlot(STobject, group.by = "Cancer", pt.size.factor = 0.65, image.alpha = 1, stroke = 0.25)  + 
  scale_fill_manual(labels = c("TME-distant", "TME-proximal"), values = c("#74ADD1", "#313695"), na.translate = F) + 
  theme(legend.key = element_blank(), legend.direction = "horizontal", legend.position = "top", legend.title = element_blank(), legend.text=element_text(size=35)) + 
  guides(fill = guide_legend(override.aes = list(size = 4)))
p
pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/figures/supp/Cancer_", name, cancer_factors, ".pdf"), width = 9, height = 7)
plot(p)
dev.off()



STobject@meta.data[,paste0("H",1:13)] <- scale(STobject@meta.data[,paste0("H",1:13)])

H_modules <- STobject@meta.data[,c("Cancer", paste0("H", 1:13))]
H_modules[2:14] <- scale(H_modules[2:14])
colnames(H_modules[2:14]) <- hallmark_names
mean(H_modules$H12)
sd(H_modules$H12)

library(reshape2)
H_modules_resh <- melt(H_modules,id.vars='Cancer', measure.vars=paste0("H", 1:13))
colnames(H_modules_resh)[1:2] <- c("Cancer regions", "Hallmark")
H_modules_resh <- H_modules_resh[!is.na(H_modules_resh$`Cancer regions`),]
H_modules_resh <- H_modules_resh[H_modules_resh$Hallmark %in% paste0("H", c(2,4,8,9,11,12)),]

H_modules_resh$Hallmark <- sapply(as.character(H_modules_resh$Hallmark), function(x) {
  n <- substr(x, 2, stop = nchar(as.character(x)))
  hallmark_names[as.numeric(n)]
})

H2 <- H_modules_resh

ggboxplot(H_modules_resh, x = "Hallmark", y = "value", ylab = "Hallmark scores", merge = "flip",
          ggtheme = theme_classic(), fill = "Cancer regions", palette =  c( "#313695","#74ADD1"),
          legend = "top",
          title = "",
          font.legend = c(15, "plain", "black"),
          font.tickslab = c(15,"black"),
          font.main = c(15, "plain", "black"),
          font.x = c(15, "plain", "black"),
          font.y = c(15, "plain", "black") 
) + stat_compare_means() + theme(axis.text.x = element_text(angle = -45, hjust = 0))


for (hallmark in hallmark_names[c(2,4,8,9,11,12)]) {
  group1 <- H_modules_resh$value[H_modules_resh$`Cancer regions`=="TME-distant" & H_modules_resh $Hallmark == hallmark]
  group2 <- H_modules_resh$value[H_modules_resh$`Cancer regions`=="TME-proximal" & H_modules_resh $Hallmark == hallmark]
  t <- wilcox.test(group1, 
                   group2, 
                   alternative = "two.sided")
  diff_means <- mean(group2) - mean(group1)
  row <- c(name, cancer_factors,  hallmark, diff_means, t$p.value)
  summary_comparisions <- rbind(summary_comparisions, row)
}


#TME
min_dist <- apply(distances,2,min)

inner <- names(min_dist[min_dist >= quantile(min_dist, 0.75)])
outer <- names(min_dist[min_dist <= quantile(min_dist, 0.25)])


STobject$TME <- NA
STobject$TME[inner] <- "Cancer-distant"
STobject$TME[outer] <- "Cancer-proximal"

p <- SpatialDimPlot(STobject, group.by = "TME", pt.size.factor = 0.65, image.alpha = 1)  + 
  scale_fill_manual(labels = c("Cancer-distant", "Cancer-proximal"), values = c("#F46D43", "#A50026"), na.translate = F) + 
  theme(legend.key = element_blank(), legend.direction = "horizontal", legend.position = "top", legend.title = element_blank(), legend.text=element_text(size=35)) + 
  guides(fill = guide_legend(override.aes = list(size = 4)))
p
pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/figures/supp/TME_", name, cancer_factors, ".pdf"), width = 9, height = 7)
plot(p)
dev.off()


H_modules <- STobject@meta.data[,c("TME", paste0("H", 1:13))]
H_modules[2:14] <- scale(H_modules[2:14])
mean(H_modules$H12)
sd(H_modules$H12)

H_modules_resh <- melt(H_modules,id.vars='TME', measure.vars=paste0("H", 1:13))
colnames(H_modules_resh)[1:2] <- c("TME regions", "Hallmark")
H_modules_resh <- H_modules_resh[!is.na(H_modules_resh$`TME regions`),]
H_modules_resh <- H_modules_resh[H_modules_resh$Hallmark %in% paste0("H", c(1,3,5,6,7,10,13)),]
H_modules_resh$Hallmark <- sapply(as.character(H_modules_resh$Hallmark), function(x) {
  n <- substr(x, 2, stop = nchar(as.character(x)))
  hallmark_names[as.numeric(n)]
})

H1 <- H_modules_resh

ggboxplot(H_modules_resh, x = "TME regions", y = "value", ylab = "Hallmark scores", merge = "flip",
          ggtheme = theme_classic(), fill = "TME regions", palette = c("#F46D43", "#A50026"),
          legend = "top",
          title = "",
          font.legend = c(15, "plain", "black"),
          font.tickslab = c(15,"black"),
          font.main = c(15, "plain", "black"),
          font.x = c(15, "plain", "black"),
          font.y = c(15, "plain", "black"), facet.by = "Hallmark"
) + stat_compare_means()

for (hallmark in hallmark_names[c(1,3,5,6,7,10,13)]) {
  group1 <- H_modules_resh$value[H_modules_resh$`TME regions`=="Cancer-distant" & H_modules_resh $Hallmark == hallmark]
  group2 <- H_modules_resh$value[H_modules_resh$`TME regions`=="Cancer-proximal" & H_modules_resh $Hallmark == hallmark]
  t <- wilcox.test(group1, 
                   group2, 
                   alternative = "two.sided")
  diff_means <- mean(group2) - mean(group1)
  row <- c(name, cancer_factors,  hallmark, diff_means, t$p.value)
  summary_comparisions <- rbind(summary_comparisions, row)
}


#summary_comparisions <- read.table("Desktop/comp.txt", header = T, sep = "\t")
#write.table(summary_comparisions, "Desktop/comp.txt", sep = "\t")

colnames(H1)[1] <- "Tumor regions"
colnames(H2)[1] <- "Tumor regions"
H <- rbind(H1, H2)
H

p <- ggboxplot(H, x = "Hallmark", y = "value", ylab = "Scaled Hallmark scores", merge = "flip",
          ggtheme = theme_classic(), fill = "Tumor regions", palette = c("#F46D43", "#A50026","#74ADD1", "#313695"),
          legend = "top",
          title = "",
          font.legend = c(20, "plain", "black"),
          font.tickslab = c(15,"black"),
          font.main = c(15, "plain", "black"),
          font.x = c(15, "plain", "black"),
          font.y = c(15, "plain", "black"),font.label = c(35, "plain", "black"),
)  + theme(axis.text.x = element_text(angle = -45, hjust = 0, size = 15), plot.margin = margin(0.5,6.5,0.5,0.5, "cm"))  

p
pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/figures/supp/boxplot_", name, cancer_factors,".pdf"), width = 13, height = 8)
plot(p)
dev.off()
