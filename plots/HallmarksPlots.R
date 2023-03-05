##################################################
## Project: Cancer Hallmarks
## Script purpose: Reproduce plots
## Date: 28/12/2022
## Author: Sergi Cervilla* & Mustafa Sibai*
##################################################
# Plot outline

## Highlights
### P1
#### P1c: Number of pathways per Hallmark (barplot)
#### P1d: Number of genes per Hallmark (barplot)
#### P1e: Number of genes shared by Hallmarks (barplot)
#### P1f: Hallmark activity and marker expression plots
##### Pancreas
##### Liver
##### Kidney
##### Colorectal
##### Breast
##### Brain
##### Bladder
### P2 
#### P2a2: ESTIMATE score in colorectal tissue 
#### P2a3: ESTIMATE clusters in colorectal tissue 
#### P2a4: ESTIMATE validation boxplots
#### P2b: Pan-Cancer heatmap
### P3
#### P3b1: ESTIMATEvsNeighborhood (scatterplot)
#### P3b2: Residual values in colorectal tissue
#### P3c: R-squared of ESTIMATE+Neighborhood, all hallmarks (barplot)
#### P3d: R-squared of ESTIMATE+Neighborhood, split by hallmark (circular barplot)
#### P3e: R-squared of ESTIMATE+Neighborhood, split by hallmark and compartment (layered circular barplot)
### P5
#### P5a: Random Forest for Cancer Hallmarks results (circos plot)
#### P5b: Random Forest for Cancer Hallmarks contributions and dependencies (boxplot)
#### P5l: Cluster plots for SCD legend
### P6
#### P6a: Random Forest for TME Hallmarks results (circos plot)
#### P6b: Random Forest for TME Hallmarks contributions and dependencies (boxplot)
#### P6c2: M3 clusters
#### P6c3: M3 H3-H11
#### P6c4 : Random H3-H10
### P7
#### P7a: Clonal hallmark activity difference (heatmap)
#### P7b1: Clonal clusters 
#### P7b2: Hallmark activity across clones (boxplot)
#### P7c: Pathway activity plots in M3





# Code
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
#path to SamplesMetadata.R file
source("../utils/SamplesMetadata.R")
#path to PlottingMod.R file
source("../utils/PlottingMod.R")

## Highlights
### P1
#### P1c (GeneCollection.R)
paths_freq <- data.frame(table(paths_final$hallmark))
colnames(paths_freq)[1:2] <- c("Hallmark", "Freq")

paths_freq$Hallmark <- factor(paths_freq$Hallmark, levels = paste0("H", 1:13))

# in order to put Hallmark names, arrange by Hallmark order
paths_freq <- dplyr::arrange(paths_freq, Hallmark)
paths_freq$Hallmark <- hallmark_names
# arrange based on frequency
paths_freq <- dplyr::arrange(paths_freq, desc(Freq))

ggbarplot(paths_freq, x = "Hallmark", y = "Freq", fill = "black", color = "white", width = 0.9) + rotate_x_text(angle = -45, vjust = 2, hjust = 0) + ylab("n.Pathways")

#### P1d (GeneCollection.R)
H.genes.tbl <- data.frame(n_genes = table(H.genes$H))
colnames(H.genes.tbl)[1:2] <- c("Hallmark", "n")

H.genes.tbl$Hallmark <- factor(H.genes.tbl$Hallmark, levels = paste0("H", 1:13))
H.genes.tbl <- dplyr::arrange(H.genes.tbl, Hallmark)
H.genes.tbl$Hallmark <- hallmark_names
# arrange based on pathway frequency
H.genes.tbl$Hallmark <- factor(H.genes.tbl$Hallmark, levels = paths_freq$Hallmark)

ggbarplot(H.genes.tbl, x = "Hallmark", y = "n", fill = "black", color = "white", width = 0.9) + rotate_x_text(angle = -45, vjust = 2, hjust = 0) +  ylab("n.Genes") 


#### P1e (GeneCollection.R) 
H.genes.tbl <- data.frame(n_genes = table(H.genes$gene))
colnames(H.genes.tbl)[1:2] <- c("genes", "HAGs")

H.genes.tbl <- data.frame(table(H.genes.tbl$HAGs))
colnames(H.genes.tbl)[1:2] <- c("HAGs", "Freq")

ggbarplot(H.genes.tbl, x = "HAGs", y = "Freq", fill = "black", color = "white", width = 0.9) + ylab("n.Genes") + xlab("n. associated Hallmarks per Gene")

#### P1f
#load spot_info (ComputeNeighborsScores.R)
##### Pancreas
sample <- "P306"
#load enhanced object with imputed genes
sce <- readRDS("")
sce <- sce[,-1]
h <- "H11"
colData(sce)[,h] <- spot_info[spot_info$sample==sample, h]

v <- .make_triangle_subspots(colData(sce), fill = h)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(hallmark_names_list[[h]]) + theme(plot.title = element_text(size=35, hjust = 0.5),
                                                                                                   legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
#save file 
pdf(, width = 10, height = 10) 
plot(p)
dev.off()
g <- "TXN"
sce$gene <- assay(sce)[g,]
v <- .make_triangle_subspots(colData(sce), fill = gene)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(g) + theme(plot.title = element_text(size=50, hjust = 0.5),
                                                                            legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
#save file 
pdf("", width = 10, height = 10) 
plot(p)
dev.off()

##### Liver
sample <- "HCC2T"
#load enhanced object with imputed genes
sce <- readRDS("")
sce <- sce[,-1]
h <- "H2"
colData(sce)[,h] <- spot_info[spot_info$sample==sample, h]

v <- .make_triangle_subspots(colData(sce), fill = h)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(hallmark_names_list[[h]]) + theme(plot.title = element_text(size=35, hjust = 0.5),
                                                                                                   legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
#save file 
pdf("", width = 10, height = 10) 
plot(p)
dev.off()
g <- "CDK1"
sce$gene <- assay(sce)[g,]
v <- .make_triangle_subspots(colData(sce), fill = gene)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(g) + theme(plot.title = element_text(size=50, hjust = 0.5),
                                                                            legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
#save file 
pdf("", width = 10, height = 10) 
plot(p)
dev.off()

##### Kidney
sample <- "C21"
#load enhanced object with imputed genes
sce <- readRDS("")
sce <- sce[,-1]
h <- "H1"
colData(sce)[,h] <- spot_info[spot_info$sample==sample, h]

v <- .make_triangle_subspots(colData(sce), fill = h)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(hallmark_names_list[[h]]) + theme(plot.title = element_text(size=35, hjust = 0.5),
                                                                                                   legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
#save file 
pdf("", width = 10, height = 10) 
plot(p)
dev.off()
g <- "AKT1"
sce$gene <- assay(sce)[g,]
v <- .make_triangle_subspots(colData(sce), fill = gene)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(g) + theme(plot.title = element_text(size=50, hjust = 0.5),
                                                                            legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
#save file 
pdf("", width = 10, height = 10) 
plot(p)
dev.off()

##### Colorectal
sample <- "Intestine"
#load enhanced object with imputed genes
sce <- readRDS("")
sce <- sce[,-1]
h <- "H8"
colData(sce)[,h] <- spot_info[spot_info$sample==sample, h]

v <- .make_triangle_subspots(colData(sce), fill = h)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(hallmark_names_list[[h]]) + theme(plot.title = element_text(size=35, hjust = 0.5),
                                                                                                   legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
#save file 
pdf("", width = 10, height = 10) 
plot(p)
dev.off()
g <- "BRCA1"
sce$gene <- assay(sce)[g,]
v <- .make_triangle_subspots(colData(sce), fill = gene)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(g) + theme(plot.title = element_text(size=50, hjust = 0.5),
                                                                            legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
#save file 
pdf("", width = 10, height = 10) 
plot(p)
dev.off()

##### Breast
sample <- "DuctalFFPE"
#load enhanced object with imputed genes
sce <- readRDS("")
sce <- sce[,-1]
h <- "H12"
colData(sce)[,h] <- spot_info[spot_info$sample==sample, h]

v <- .make_triangle_subspots(colData(sce), fill = h)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(hallmark_names_list[[h]]) + theme(plot.title = element_text(size=35, hjust = 0.5),
                                                                                                   legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
#save file 
pdf("", width = 10, height = 10) 
plot(p)
dev.off()
g <- "EZH2"
sce$gene <- assay(sce)[g,]
v <- .make_triangle_subspots(colData(sce), fill = gene)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(g) + theme(plot.title = element_text(size=50, hjust = 0.5),
                                                                            legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
#save file 
pdf("", width = 10, height = 10) 
plot(p)
dev.off()
##### Glioblastoma
sample <- "Glioblastoma"
#load enhanced object with imputed genes
sce <- readRDS("")
sce <- sce[,-1]
h <- "H10"
colData(sce)[,h] <- spot_info[spot_info$sample==sample, h]

v <- .make_triangle_subspots(colData(sce), fill = h)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(hallmark_names_list[[h]]) + theme(plot.title = element_text(size=35, hjust = 0.5),
                                                                                                   legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
#save file 
pdf("", width = 10, height = 10) 
plot(p)
dev.off()
g <- "VDAC1"
sce$gene <- assay(sce)[g,]
v <- .make_triangle_subspots(colData(sce), fill = gene)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(g) + theme(plot.title = element_text(size=50, hjust = 0.5),
                                                                            legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
#save file 
pdf("", width = 10, height = 10) 
plot(p)
dev.off()
##### Bladder
sample <- "DU2"
#load enhanced object with imputed genes
sce <- readRDS("")
sce <- sce[,-1]
h <- "H6"
colData(sce)[,h] <- spot_info[spot_info$sample==sample, h]

v <- .make_triangle_subspots(colData(sce), fill = h)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(hallmark_names_list[[h]]) + theme(plot.title = element_text(size=35, hjust = 0.5),
                                                                                                   legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
#save file 
pdf("", width = 10, height = 10) 
plot(p)
dev.off()
g <- "VCAN"
sce$gene <- assay(sce)[g,]
v <- .make_triangle_subspots(colData(sce), fill = gene)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(g) + theme(plot.title = element_text(size=50, hjust = 0.5),
                                                                            legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
#save file 
pdf("", width = 10, height = 10) 
plot(p)
dev.off()

### P2 
#image
hires <- readRDS("../image/Colorectal.rds")
#spot location
v <- readRDS("../bayes/Colorectal.rds")
v <- v[v$spot != "subspot_1.1",]
#sub-spot data
sub_data <- read.table("../output_df/Colorectal.txt")
cluster <- c("5" = "lightgoldenrod1", "4" = "lightgoldenrod3", "3" = "lightpink2", "2"= "orchid3", "1"= "orchid4")

#### P2a2
v$fill <- sub_data[v$spot,"estimate"]
#save file 
pdf("", width = 10, height = 10) 
hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_gradientn(colours = rev(cluster), breaks=c(min(v$fill), max(v$fill)),labels=c("Cancer pure","TME pure")) + labs(fill="") + 
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.key.width = unit(118, "point"),
        legend.text = element_text(size = 20), plot.title = element_text(size = 45, hjust = 0.5)) +
  ggtitle("ESTIMATE score")
dev.off()

#### P2a3
v$fill <- sub_data[v$spot,"clusters"]
#save file 
pdf("", width = 10, height = 10) 
hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~factor(fill))) +  theme_void() + coord_equal() +
  scale_fill_manual(values = rev(c("lightgoldenrod1", "lightgoldenrod3", "lightpink2",  "orchid3",  "orchid4"))) + labs(fill="") + 
  ggtitle("ESTIMATE clusters") + theme(plot.title = element_text(size = 45, hjust = 0.5), legend.key.size = unit(30,"points"),
                                       legend.text = element_text(size = 20))
dev.off()

#### P2a4
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
#save file 
pdf("", width = 12, height = 8) 
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
#save file 
pdf("", width = 12, height = 8) 
p1/p2 
dev.off()

#### P2b
# (PanCancerHeatmap.R)

### P3
#### P3b1
d <- data.frame(estimate=sub_data$estimate, score=sub_data$score)
model <- lm(estimate~score, data = d)
sub_data$residuals <- model$residuals
#save file 
pdf("", width = 10, height = 7) 
ggplot(sub_data, aes(x=score, y=estimate, col=residuals)) + geom_point() + geom_smooth(method = "lm", color="black") +
  scale_color_gradient2(low = "navy", mid = "gainsboro", high = "firebrick3") + theme_classic() + labs(x="NEIGHBORHOOD score", y="ESTIMATE score", col="Residuals") + 
  scale_x_continuous(labels = c("Cancer", "TME"), breaks = c(-50, 100)) + 
  scale_y_continuous(labels = c("Cancer", "TME"), breaks = c(-5000, 10000)) + 
  theme(axis.ticks = element_blank(), axis.text = element_text(size = 20), axis.title = element_text(size = 25), legend.position = "none") 
dev.off()

#### P3b2
v$fill <- sub_data[v$spot,"residuals"]
#save file 
pdf("", width = 10, height = 10) 
hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_gradient2(low = "navy", mid = "gainsboro", high = "firebrick3") + labs(fill="") + ggtitle("Residuals") + 
  theme(plot.title = element_text(size = 35, hjust = 0.5), legend.key.height = unit(115, "points"), 
        legend.key.width = unit(25, "points"), legend.text = element_text(size=20))
dev.off()

#### P3c
#### P3d
#### P3e
### P5
#### P5a
# (CancerCircos.R)

#### P5b
# (CancerCircos.R)

#### P5l: Cluster plots for SCD legend
sample <- "M1"
#single cell experiment enhanced object without imputed genes
sce <- readRDS("")
sce <- sce[,-1]
sce$clusters <- spot_info[spot_info$sample==sample, "clusters"]
sce$compartments <- cut(sce$clusters, breaks = c(0,2.5,3.5,5), labels = c("Cancer", "Buffer", "TME"))
v <- .make_triangle_subspots(colData(sce), fill = "compartments")
#spot location
ref_v <- readRDS("")
#image
hires <- readRDS("")
for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}
#save file 
pdf("", width = 7, height = 7) 
hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~as.factor(fill))) +  theme_void() +
  labs(fill="Compartments") +  scale_fill_manual(values = rev(c("lightgoldenrod1",  "lightpink2",  "orchid3")))  + theme(legend.position = "bottom")
dev.off()

sample <- "DU13"
#single cell experiment enhanced object without imputed genes
sce <- readRDS("")
sce <- sce[,-1]
sce$clusters <- spot_info[spot_info$sample==sample, "clusters"]
sce$compartments <- cut(sce$clusters, breaks = c(0,2.5,3.5,5), labels = c("Cancer", "Buffer", "TME"))
v <- .make_triangle_subspots(colData(sce), fill = "compartments")
#spot location
ref_v <- readRDS("")
#image
hires <- readRDS("")
for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}
#save file 
pdf("", width = 7, height = 7) 
hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~as.factor(fill))) +  theme_void() +
  labs(fill="Compartments") +  scale_fill_manual(values = rev(c("lightgoldenrod1",  "lightpink2",  "orchid3")))  + theme(legend.position = "bottom")
dev.off()


### P6
#### P6a
# (TMECircos.R)
#### P6b
# (TMECircos.R)
#### P6c2: M3 clusters

#### P6c3
#load sub-spot data of all samples
spot_info <- read.table("")
sample <- "M3"
print(sample)
h_cancer <- "H11"
h_tme <- "H3"
#single cell experiment enhanced object without imputed genes
sce <- readRDS("")
sce <- sce[,-1]
colData(sce)[, c(h_cancer, h_tme, "clusters")] <- spot_info[spot_info$sample==sample, c(h_cancer, h_tme, "clusters")]

v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill = h_cancer)
v2 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(4,5),], fill = h_tme)
v3 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(3),], fill = "clusters")

#spot location
ref_v <- readRDS("")
#image
hires <- readRDS("")
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
#save file 
pdf("", width = 7, height = 7) 
hires + geom_polygon(data=v,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn(h_cancer, colours = viridisLite::rocket(1000, alpha = 1, begin = 0, end = 1, direction = 1)[200:1000]) + 
  ggnewscale::new_scale("fill") + geom_polygon(data=v2,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn(h_tme, colours = viridisLite::mako(1000, alpha = 1, begin = 0, end = 1, direction = 1)[200:1000]) + 
  ggnewscale::new_scale("fill") + geom_polygon(data=v3,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~as.factor(fill))) +
  scale_fill_manual("Buffer", values = "gray41") 
dev.off()
#### P6c4 
#load sub-spot data of all samples
spot_info <- read.table("Desktop/df_cnv.txt")
sample <- "M3"
print(sample)
h_cancer <- "H11"
h_tme <- "H3"
#single cell experiment enhanced object without imputed genes
sce <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/M3_sce_enhanced.rds")
sce <- sce[,-1]
colData(sce)[, c(h_cancer, h_tme, "clusters")] <- spot_info[spot_info$sample==sample, c(h_cancer, h_tme, "clusters")]

v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2, 3),], fill = h_cancer)
v2 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(4,5),], fill = h_tme)

#spot location
ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/M3.rds")
#image
hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/M3.rds")
for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}



#save file 
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/Final/M3.pdf", width = 7, height = 7) 
hires + geom_polygon(data=v,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn(h_cancer, colours = viridisLite::rocket(1000, alpha = 1, begin = 0, end = 1, direction = 1)[200:1000]) + 
  ggnewscale::new_scale("fill") + geom_polygon(data=v2,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn(h_tme, colours = viridisLite::mako(1000, alpha = 1, begin = 0, end = 1, direction = 1)[200:1000]) 
dev.off()
### P7
#### P7a 
#load df_diff (CNV_experiment.R)
col_fun = colorRamp2(c(0, max(df_diff)), c("white", "red"))
tissue_annot <- sapply(rownames(df_diff), function(x) {annotation_tissue[[x]]})
col_tissue_map <- setNames(Seurat::DiscretePalette(15, palette = "polychrome") [c(1:3, 5, 6:11)], sort(unique(tissue_annot)))

row_ha = HeatmapAnnotation(`Tumor type` = tissue_annot,
                           col=list(`Tumor type` = col_tissue_map))
palette <- do.call(c,color_codes[c(2,4,8,10,11,12)])
names(palette) <- NULL
top_ha = rowAnnotation(" " = anno_boxplot(df_diff, height = unit(2, "cm"),gp = gpar(fill = palette)))
#save file 
pdf("", width = 14, height = 3) 
Heatmap(t(df_diff), name = "Difference", border = "black", col=col_fun, show_row_names =T, 
        show_column_names =F,
        top_annotation =  row_ha,
        left_annotation = top_ha, 
        rect_gp = gpar(col = "gray", lwd = 0.1))
dev.off()
#### P7b1 
samples <- c("M3", "Intestine")
cnv_clusters <- list(M3=1:2, Intestien=0:2, HCC2T=0:2, OV4A=0)
for (sample in samples) {
  #single cell experiment enhanced object without imputed genes
  sce <- readRDS("")
  sce <- sce[,-1]
  sce$clusters <- spot_info[spot_info$sample==sample, "clusters"]
  sce$cnv <- cnv[cnv$sample == sample, "cnv_cluster"]
  v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2) & sce$cnv %in% cnv_clusters[[sample]],], fill = "cnv")
  #spot location
  ref_v <- readRDS("")
  #image
  hires <- readRDS("")
  for (spot in unique(v$spot)) {
    v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
    v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
  }
  #save file 
  pdf("", width = 7, height = 7) 
  hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~as.factor(fill))) +  theme_void() +
    labs(fill="CNV clones") +  scale_fill_manual(values = c("#E69F00", "#0072B2", "#009E73"))
  dev.off()
}
#### P7b2
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
  #save file 
  pdf("", width = 14, height = 7) 
  grid::grid.draw(gp)
  dev.off()
  
}
#### P7c 
#load paths from PathScore (CNV_experiment.R)
paths <- read.table("Downloads/M3.txt", sep = "\t", check.names = F)
paths[, 15:ncol(paths)] <- scale(paths[, 15:ncol(paths)])
#load sub-spot data of all samples (same as spot_info)
cnv <- read.table("")
cnv <- cnv %>% filter(sample == "M3")
paths[cnv$spotid, "cnv_cluster"] <- cnv$cnv_cluster


pathways <- c("H10_Glycolysis", "H11_Senescence-Associated Secretory Phenotype (SASP)",
              "H2_p53 pathway", "H4_Packaging Of Telomere Ends", "H8_G2/M DNA damage checkpoint", "H12_HATs acetylate histones")


sample <- "M3"
#single cell experiment enhanced object without imputed genes
sce <- readRDS("")
sce <- sce[,-1]
pathway <- pathways[2]
colData(sce)[, c(pathways, "clusters", "cnv_cluster")] <- paths[, c(pathways, "estimate.cluster", "cnv_cluster")]
sce$empty <- "empty"
#spot location
ref_v <- readRDS("")
#image
hires <- readRDS("")
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
  #save file 
  pdf("", width = 7, height = 7) 
  hires + geom_polygon(data=v,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
    scale_fill_gradientn("Pathway activity", colours = viridisLite::rocket(1000, alpha = 1, begin = 0, end = 1, direction = 1)[200:1000]) + 
    ggnewscale::new_scale("fill") + geom_polygon(data=v3,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~as.factor(fill))) +
    scale_fill_manual("Buffer", values = "gray41") + theme(plot.title = element_text(hjust = 0.5, size=35),
                                                           legend.position = "none") +  ggtitle("Senescence-Associated\n Secretory Phenotype (SASP)")
  ggtitle(str_split(pathway, pattern = "_", simplify = T)[1,2])
  dev.off()
}




###########

#load sub-spot data of all samples
spot_info <- read.table("Desktop/df_cnv.txt")
sample <- "Co1"
h_cancer <- "H4"
#single cell experiment enhanced object without imputed genes
sce <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/Co1_sce_enhanced.rds")
sce <- sce[,-1]
colData(sce)[, c(h_cancer, "clusters")] <- spot_info[spot_info$sample==sample, c(h_cancer, "clusters")]

v <- .make_triangle_subspots(colData(sce), fill = h_cancer)

#spot location
ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/Co1.rds")
#image
hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/Co1.rds")
for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}


#save file 
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/Final/Fig2_H4.pdf", width = 7, height = 7) 
hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn(h_cancer, colours = viridisLite::rocket(1000, alpha = 1, begin = 0, end = 1, direction = 1)[200:1000])
dev.off()

spot_info <- read.table("Desktop/df_cnv.txt")
sample <- "Co1"
h_cancer <- "H13"
#single cell experiment enhanced object without imputed genes
sce <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/Co1_sce_enhanced.rds")
sce <- sce[,-1]
colData(sce)[, c(h_cancer, "clusters")] <- spot_info[spot_info$sample==sample, c(h_cancer, "clusters")]

v <- .make_triangle_subspots(colData(sce), fill = h_cancer)

#spot location
ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/Co1.rds")
#image
hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/Co1.rds")
for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}


#save file 
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/Final/Fig2_H13.pdf", width = 7, height = 7) 
hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn(h_cancer, colours = viridisLite::rocket(1000, alpha = 1, begin = 0, end = 1, direction = 1)[200:1000])
dev.off()


##HALLMARKS 
