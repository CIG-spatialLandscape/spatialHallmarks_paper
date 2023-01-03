library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

source("Desktop/IJC/scripts/HallmarksPaper/utils/SamplesMetadata.R")
source("Desktop/IJC/datasets/IGTP/figuresPaper/scripts/utilities/BayesSpace_functions.r")
### Highlight 1 ###
#Pancreas
sample <- "P306"
sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
sce <- sce[,-1]
h <- "H11"
colData(sce)[,h] <- spot_info[spot_info$sample==sample, h]

v <- .make_triangle_subspots(colData(sce), fill = h)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(hallmark_names_list[[h]]) + theme(plot.title = element_text(size=35, hjust = 0.5),
                                                                                                   legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H1/", sample, "_", h, ".pdf"), width = 10, height = 10) 
plot(p)
dev.off()
g <- "TXN"
sce$gene <- assay(sce)[g,]
v <- .make_triangle_subspots(colData(sce), fill = gene)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(g) + theme(plot.title = element_text(size=50, hjust = 0.5),
                                                                                                   legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H1/", sample, "_", g, ".pdf"), width = 10, height = 10) 
plot(p)
dev.off()
#Liver
sample <- "HCC2T"
sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
sce <- sce[,-1]
h <- "H2"
colData(sce)[,h] <- spot_info[spot_info$sample==sample, h]

v <- .make_triangle_subspots(colData(sce), fill = h)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(hallmark_names_list[[h]]) + theme(plot.title = element_text(size=35, hjust = 0.5),
                                                                                                   legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H1/", sample, "_", h, ".pdf"), width = 10, height = 10) 
plot(p)
dev.off()
g <- "CDK1"
sce$gene <- assay(sce)[g,]
v <- .make_triangle_subspots(colData(sce), fill = gene)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(g) + theme(plot.title = element_text(size=50, hjust = 0.5),
                                                                            legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H1/", sample, "_", g, ".pdf"), width = 10, height = 10) 
plot(p)
dev.off()
#Kidney
sample <- "C21"
sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
sce <- sce[,-1]
h <- "H1"
colData(sce)[,h] <- spot_info[spot_info$sample==sample, h]

v <- .make_triangle_subspots(colData(sce), fill = h)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(hallmark_names_list[[h]]) + theme(plot.title = element_text(size=35, hjust = 0.5),
                                                                                                   legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H1/", sample, "_", h, ".pdf"), width = 10, height = 10) 
plot(p)
dev.off()
g <- "AKT1"
sce$gene <- assay(sce)[g,]
v <- .make_triangle_subspots(colData(sce), fill = gene)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(g) + theme(plot.title = element_text(size=50, hjust = 0.5),
                                                                            legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H1/", sample, "_", g, ".pdf"), width = 10, height = 10) 
plot(p)
dev.off()
#Colorectal
sample <- "Intestine"
sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
sce <- sce[,-1]
h <- "H8"
colData(sce)[,h] <- spot_info[spot_info$sample==sample, h]

v <- .make_triangle_subspots(colData(sce), fill = h)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(hallmark_names_list[[h]]) + theme(plot.title = element_text(size=35, hjust = 0.5),
                                                                                                   legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H1/", sample, "_", h, ".pdf"), width = 10, height = 10) 
plot(p)
dev.off()
g <- "BRCA1"
sce$gene <- assay(sce)[g,]
v <- .make_triangle_subspots(colData(sce), fill = gene)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(g) + theme(plot.title = element_text(size=50, hjust = 0.5),
                                                                            legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H1/", sample, "_", g, ".pdf"), width = 10, height = 10) 
plot(p)
dev.off()
#Breast
sample <- "DuctalFFPE"
sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
sce <- sce[,-1]
h <- "H12"
colData(sce)[,h] <- spot_info[spot_info$sample==sample, h]

v <- .make_triangle_subspots(colData(sce), fill = h)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(hallmark_names_list[[h]]) + theme(plot.title = element_text(size=35, hjust = 0.5),
                                                                                                   legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H1/", sample, "_", h, ".pdf"), width = 10, height = 10) 
plot(p)
dev.off()
g <- "EZH2"
sce$gene <- assay(sce)[g,]
v <- .make_triangle_subspots(colData(sce), fill = gene)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(g) + theme(plot.title = element_text(size=50, hjust = 0.5),
                                                                            legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H1/", sample, "_", g, ".pdf"), width = 10, height = 10) 
plot(p)
dev.off()
#Glioblastoma
sample <- "Glioblastoma"
sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
sce <- sce[,-1]
h <- "H10"
colData(sce)[,h] <- spot_info[spot_info$sample==sample, h]

v <- .make_triangle_subspots(colData(sce), fill = h)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(hallmark_names_list[[h]]) + theme(plot.title = element_text(size=35, hjust = 0.5),
                                                                                                   legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H1/", sample, "_", h, ".pdf"), width = 10, height = 10) 
plot(p)
dev.off()
g <- "VDAC1"
sce$gene <- assay(sce)[g,]
v <- .make_triangle_subspots(colData(sce), fill = gene)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(g) + theme(plot.title = element_text(size=50, hjust = 0.5),
                                                                            legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H1/", sample, "_", g, ".pdf"), width = 10, height = 10) 
plot(p)
dev.off()
#Bladder
sample <- "DU2"
sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
sce <- sce[,-1]
h <- "H6"
colData(sce)[,h] <- spot_info[spot_info$sample==sample, h]

v <- .make_triangle_subspots(colData(sce), fill = h)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(hallmark_names_list[[h]]) + theme(plot.title = element_text(size=35, hjust = 0.5),
                                                                                                   legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H1/", sample, "_", h, ".pdf"), width = 10, height = 10) 
plot(p)
dev.off()
g <- "VCAN"
sce$gene <- assay(sce)[g,]
v <- .make_triangle_subspots(colData(sce), fill = gene)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(g) + theme(plot.title = element_text(size=50, hjust = 0.5),
                                                                            legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H1/", sample, "_", g, ".pdf"), width = 10, height = 10) 
plot(p)
dev.off()
### Highlight 2 ###
hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/Colorectal.rds")
v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/Colorectal.rds")
v <- v[v$spot != "subspot_1.1",]
sub_data <- read.table("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/output_df/Colorectal.txt")
cluster <- c("5" = "lightgoldenrod1", "4" = "lightgoldenrod3", "3" = "lightpink2", "2"= "orchid3", "1"= "orchid4")
#continous estimate score
v$fill <- sub_data[v$spot,"estimate"]
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H2/EstimateColorectal.pdf", width = 10, height = 10)
hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_gradientn(colours = rev(cluster), breaks=c(min(v$fill), max(v$fill)),labels=c("Cancer pure","TME pure")) + labs(fill="") + 
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.key.width = unit(118, "point"),
        legend.text = element_text(size = 20), plot.title = element_text(size = 45, hjust = 0.5)) +
  ggtitle("ESTIMATE score")
dev.off()

#estimate clusters
v$fill <- sub_data[v$spot,"clusters"]
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H2/ClustersColorectal.pdf", width = 10, height = 10)
hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~factor(fill))) +  theme_void() + coord_equal() +
  scale_fill_manual(values = rev(c("lightgoldenrod1", "lightgoldenrod3", "lightpink2",  "orchid3",  "orchid4"))) + labs(fill="") + 
  ggtitle("ESTIMATE clusters") + theme(plot.title = element_text(size = 45, hjust = 0.5), legend.key.size = unit(30,"points"),
                                       legend.text = element_text(size = 20))
dev.off()
#boxplot validation
df <- readxl::read_xlsx("Downloads/Rstudiotumor(1).xlsx")
df$Others[20] <- "100"
colnames(df)[1] <- "clusters"
df_filtered <- df %>% filter(!is.na(`Tumor %`)) %>% select(clusters,`Tumor %`, `Stroma %`, `Necrosis %`)
p1 <- ggplot(df_filtered, aes(x=clusters, y=`Tumor %`, fill=clusters)) + geom_boxplot() +
  scale_fill_manual(values = c("cluster 5" = "lightgoldenrod1", "cluster 4" = "lightgoldenrod3", 
                               "cluster 3" = "lightpink2", "cluster 2"= "orchid3", "cluster 1"= "orchid4")) + 
  theme_pubr(base_size = 20) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 35)) + 
  labs(x="", y="Cancer cells %") + scale_x_discrete(labels=1:5) + ggtitle("Manual annotation of 17 samples")
p2 <- ggplot(df_filtered, aes(x=clusters, y=`Stroma %`, fill=clusters)) + geom_boxplot() +
  scale_fill_manual(values = c("cluster 5" = "lightgoldenrod1", "cluster 4" = "lightgoldenrod3", 
                               "cluster 3" = "lightpink2", "cluster 2"= "orchid3", "cluster 1"= "orchid4")) + 
  theme_pubr(base_size = 20) + theme(legend.position = "none") + labs(x="", y="Stroma cells %") + scale_x_discrete(labels=1:5)
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H2/BoxplotValidation.pdf", width = 12, height = 8)
p1/p2 
dev.off()
p1 <- ggplot(df_filtered, aes(x=clusters, y=`Tumor %`, fill=clusters)) + geom_boxplot(outlier.shape = 1) +
  scale_fill_manual(values = c("cluster 5" = "lightgoldenrod1", "cluster 4" = "lightgoldenrod3", 
                               "cluster 3" = "lightpink2", "cluster 2"= "orchid3", "cluster 1"= "orchid4")) + 
  theme_pubr(base_size = 20) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 35)) + 
  labs(x="", y="Cancer cells %") + scale_x_discrete(labels=1:5) + ggtitle("Manual annotation of 17 samples") + geom_jitter(shape=20, width = 0.3)
p2 <- ggplot(df_filtered, aes(x=clusters, y=`Stroma %`, fill=clusters)) + geom_boxplot(outlier.shape = 1) +
  scale_fill_manual(values = c("cluster 5" = "lightgoldenrod1", "cluster 4" = "lightgoldenrod3", 
                               "cluster 3" = "lightpink2", "cluster 2"= "orchid3", "cluster 1"= "orchid4")) + 
  theme_pubr(base_size = 20) + theme(legend.position = "none") + labs(x="", y="Stroma cells %") + scale_x_discrete(labels=1:5) + geom_jitter(shape=20, width = 0.3)
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H2/BoxplotValidation_dots.pdf", width = 12, height = 8)
p1/p2 
dev.off()

### Highlight 3 ###
d <- data.frame(estimate=sub_data$estimate, score=sub_data$score)
model <- lm(estimate~score, data = d)
sub_data$residuals <- model$residuals

#estimate vs neighborhood
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H3/Residuals.pdf", width = 10, height = 7)
ggplot(sub_data, aes(x=score, y=estimate, col=residuals)) + geom_point() + geom_smooth(method = "lm", color="black") +
  scale_color_gradient2(low = "navy", mid = "gainsboro", high = "firebrick3") + theme_classic() + labs(x="NEIGHBORHOOD score", y="ESTIMATE score", col="Residuals") + 
  scale_x_continuous(labels = c("Cancer", "TME"), breaks = c(-50, 100)) + 
  scale_y_continuous(labels = c("Cancer", "TME"), breaks = c(-5000, 10000)) + 
  theme(axis.ticks = element_blank(), axis.text = element_text(size = 20), axis.title = element_text(size = 25), legend.position = "none") 
dev.off()
#residuals on the tissue
v$fill <- sub_data[v$spot,"residuals"]
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H3/ResidualsTissue.pdf", width = 10, height = 10)
hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_gradient2(low = "navy", mid = "gainsboro", high = "firebrick3") + labs(fill="") + ggtitle("Residuals") + 
  theme(plot.title = element_text(size = 35, hjust = 0.5), legend.key.height = unit(115, "points"), 
        legend.key.width = unit(25, "points"), legend.text = element_text(size=20))
dev.off()
### Highlight 4 ###

### Highlight 5 ###

### Highlight 6 ###
spot_info <- read.table("Desktop/df_cnv.txt")
sample <- "M3"
print(sample)
h_cancer <- "H11"
h_tme <- "H3"
sce <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/enhanced/", sample, "_sce_enhanced.rds"))
sce <- sce[,-1]
colData(sce)[, c(h_cancer, h_tme, "clusters")] <- spot_info[spot_info$sample==sample, c(h_cancer, h_tme, "clusters")]

v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill = h_cancer)
v2 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(4,5),], fill = h_tme)
v3 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(3),], fill = "clusters")

ref_v <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/",sample,".rds"))
hires <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/",sample,".rds"))
for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}


for (spot in unique(v2$spot)) {
  v2$imagecol[v2$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v2$imagerow[v2$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}

for (spot in unique(v3$spot)) {
  v3$imagecol[v3$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v3$imagerow[v3$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}
pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H6/", sample, "_", h_cancer, "_", h_tme, ".pdf"), width = 7, height = 7)
hires + geom_polygon(data=v,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn(h_cancer, colours = viridisLite::rocket(1000, alpha = 1, begin = 0, end = 1, direction = 1)[200:1000]) + 
  ggnewscale::new_scale("fill") + geom_polygon(data=v2,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn(h_tme, colours = viridisLite::mako(1000, alpha = 1, begin = 0, end = 1, direction = 1)[200:1000]) + 
  ggnewscale::new_scale("fill") + geom_polygon(data=v3,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~as.factor(fill))) +
  scale_fill_manual("Buffer", values = "gray41") 
dev.off()
### Highlight 7 ###
spot_info <- read.table("Desktop/df_cnv.txt")
cnv <- spot_info[!is.na(spot_info$cnv_cluster), ]
cnv$annot <- "TME"
cnv$annot[cnv$clusters == 3] <- "Buffer"
cnv$annot[cnv$clusters %in% c(1,2)] <- "Cancer"
perc <- cnv %>% group_by(sample, cnv_cluster, annot) %>% 
  count() %>% group_by(sample, cnv_cluster) %>% mutate(percent=n/sum(n))  %>% 
  select(-n) %>% spread(annot, value = percent) %>% 
  mutate(BufferCancer = Buffer+Cancer) %>% filter(!is.na(cnv_cluster))
perc_filter <- perc[(perc$Cancer > 0.8 | (perc$BufferCancer > 0.65) & perc$Cancer > 0.35),] %>% filter(!is.na(sample))
sub_cnv <- data.frame()
for (sample_s in unique(perc_filter$sample)) {
  for (cluster in filter(perc_filter, sample == sample_s)) {
    sub_cnv <- rbind(sub_cnv, cnv[cnv$sample==sample_s & cnv$cnv_cluster==cluster & cnv$clusters %in% c(1,2),])
  }
}
sub_cnv <- sub_cnv %>% filter(!is.na(sample))
colnames(sub_cnv) <- colnames(cnv)
sub_cnv$cnv_cluster <- factor(sub_cnv$cnv_cluster)
test <- sub_cnv %>% group_by(sample, cnv_cluster) %>% summarise(across(H1:H13, mean)) %>% 
  group_by(sample) %>% summarise(minH1 = min(H1), maxH1=max(H1), diffH1=max(H1)-min(H1),
                                 minH2 = min(H2), maxH2=max(H2), diffH2=max(H2)-min(H2),
                                 minH3 = min(H3), maxH3=max(H3), diffH3=max(H3)-min(H3),
                                 minH4 = min(H4), maxH4=max(H4), diffH4=max(H4)-min(H4),
                                 minH5 = min(H5), maxH5=max(H5), diffH5=max(H5)-min(H5),
                                 minH6 = min(H6), maxH6=max(H6), diffH6=max(H6)-min(H6),
                                 minH7 = min(H7), maxH7=max(H7), diffH7=max(H7)-min(H7),
                                 minH8 = min(H8), maxH8=max(H8), diffH8=max(H8)-min(H8),
                                 minH9 = min(H9), maxH9=max(H9), diffH9=max(H9)-min(H9),
                                 minH10 = min(H10), maxH10=max(H10), diffH10=max(H10)-min(H10),
                                 minH11 = min(H11), maxH11=max(H11), diffH11=max(H11)-min(H11),
                                 minH12 = min(H12), maxH12=max(H12), diffH12=max(H12)-min(H12),
                                 minH13 = min(H13), maxH13=max(H13), diffH13=max(H13)-min(H13))

df_diff <- data.frame(row.names = test$sample, test[, paste0("diffH", c(2,4,8,10,11,12))])
#Heatmap
col_fun = colorRamp2(c(0, max(df_diff)), c("white", "red"))
tissue_annot <- sapply(rownames(df_diff), function(x) {annotation_tissue[[x]]})
col_tissue_map <- setNames(Seurat::DiscretePalette(15, palette = "polychrome") [c(1:3, 5, 6:11)], sort(unique(tissue_annot)))

row_ha = HeatmapAnnotation(`Tumor type` = tissue_annot,
                           col=list(`Tumor type` = col_tissue_map))
palette <- do.call(c,color_codes[c(2,4,8,10,11,12)])
names(palette) <- NULL
top_ha = rowAnnotation(" " = anno_boxplot(df_diff, height = unit(2, "cm"),gp = gpar(fill = palette)))
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H7/Heatmap.pdf", width = 14, height = 3)
Heatmap(t(df_diff), name = "Difference", border = "black", col=col_fun, show_row_names =T, 
        show_column_names =F,
        top_annotation =  row_ha,
        left_annotation = top_ha, 
        rect_gp = gpar(col = "gray", lwd = 0.1))
dev.off()
#Clusters
sample <- "OV4A"
sce <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/enhanced/",sample,"_sce_enhanced.rds"))
sce <- sce[,-1]
sce$clusters <- spot_info[spot_info$sample==sample, "clusters"]
sce$cnv <- cnv[cnv$sample == sample, "cnv_cluster"]
v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2) & sce$cnv == 0,], fill = "cnv")
ref_v <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/", sample, ".rds"))
hires <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/", sample, ".rds"))
for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}
pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H7/clusters_", sample, ".pdf"), width = 7, height = 7)
hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~as.factor(fill))) +  theme_void() +
  labs(fill="CNV clones") +  scale_fill_manual(values = c("#E69F00", "#0072B2", "#009E73"))
dev.off()
#Boxplots
samples <- c("M3", "OV4A", "HCC2T","Intestine")
for (sample_id in samples) {
  df <- spot_info[spot_info$sample == sample_id &  spot_info$cnv  %in%  perc_filter$cnv_cluster[perc_filter$sample == sample_id] & spot_info$clusters %in% 1:2, ] %>% 
    select(H2,H4,H8,H10,H11,H12, cnv_cluster) %>% melt(id.vars = "cnv_cluster")
  df$variable <- factor(df$variable, levels = c(paste0("H", c(11, 4, 12, 8, 2, 10))))
  
  p <- ggplot(df, aes(x=as.factor(variable), y=value, fill=variable)) + geom_boxplot(outlier.size = 0.25) + 
    geom_hline(aes(yintercept = 0), linetype="dashed", col="gray32") + 
    facet_grid(~as.factor(cnv_cluster), ) + 
    theme_classic() + scale_fill_manual(values = c("#05F3EB",  "#4969D4", "#890269", "#132892", "#701717", "#71189E")) + 
    labs(x="", y="Hallmark Activity") + theme(legend.position = "none", axis.text.x = element_blank(),
                                              axis.ticks.x = element_blank(), axis.text.y = element_text(size = 45),
                                              axis.title.y =  element_text(size = 55),
                                              strip.text.x = element_text(size = 0),
                                              panel.spacing = unit(3, "lines")) 
  # convert to grob
  gp <- ggplotGrob(p) # where p is the original ggplot object
  
  # assign the first 4 right-side facet strips with blue fill
  cnv_colors <- c("#E69F00", "#0072B2", "#009E73")
  gp$heights[6] <- unit(0.5, "cm")
  for (i in 1:(length(perc_filter$cnv_cluster[perc_filter$sample == sample_id]))) {
    grob.i <- grep("strip-t", gp$layout$name)[i]
    gp$grobs[[grob.i]]$grobs[[1]]$children[[1]]$gp$fill <- cnv_colors[i]
    #gp$grobs[[grob.i]]$grobs[[1]]$children[[1]]$height <- unit(2, "cm")
  }
  pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H7/boxplot_", sample_id, ".pdf"), width = 14, height = 7)
  grid::grid.draw(gp)
  dev.off()
  
}



#Pathways
paths <- read.table("Downloads/M3.txt", sep = "\t", check.names = F)
paths[, 15:ncol(paths)] <- scale(paths[, 15:ncol(paths)])
cnv <- read.table("Desktop/df_cnv.txt")
cnv <- cnv %>% filter(sample == "M3")
paths[cnv$spotid, "cnv_cluster"] <- cnv$cnv_cluster


pathways <- c("H10_Glycolysis", "H11_Senescence-Associated Secretory Phenotype (SASP)",
              "H2_p53 pathway", "H4_Packaging Of Telomere Ends", "H8_G2/M DNA damage checkpoint", "H12_HATs acetylate histones")


sample <- "M3"
sce <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/enhanced/", sample, "_sce_enhanced.rds"))
sce <- sce[,-1]
pathway <- pathways[2]
colData(sce)[, c(pathways, "clusters", "cnv_cluster")] <- paths[, c(pathways, "estimate.cluster", "cnv_cluster")]
sce$empty <- "empty"
ref_v <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/",sample,".rds"))
hires <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/",sample,".rds"))
v3 <- .make_triangle_subspots(colData(sce)[!(sce$clusters %in% c(1,2) & sce$cnv_cluster %in% 0:2),], fill = "empty")
for (spot in unique(v3$spot)) {
  v3$imagecol[v3$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v3$imagerow[v3$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}

for (pathway in pathways) {
  pathway <- pathways[2]
  print(pathway)
  v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2) & sce$cnv_cluster %in% 0:2,], fill = pathway)
  for (spot in unique(v$spot)) {
    v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
    v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
  }
  pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/Final/H7/", sample, "_", pathway, ".pdf"), bg = "white", width = 7, height = 7)
  
  hires + geom_polygon(data=v,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
    scale_fill_gradientn("Pathway activity", colours = viridisLite::rocket(1000, alpha = 1, begin = 0, end = 1, direction = 1)[200:1000]) + 
    ggnewscale::new_scale("fill") + geom_polygon(data=v3,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~as.factor(fill))) +
    scale_fill_manual("Buffer", values = "gray41") + theme(plot.title = element_text(hjust = 0.5, size=35),
                                                           legend.position = "none") +  ggtitle("Senescence-Associated\n Secretory Phenotype (SASP)")
    ggtitle(str_split(pathway, pattern = "_", simplify = T)[1,2])
  dev.off()
}
