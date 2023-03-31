##################################################
## Project: Cancer Hallmarks
## Script purpose: Reproduce main plots
## Author: Sergi Cervilla* & Mustafa Sibai*
##################################################

# Main figure plots
### P1
#### P1c: Number of genes per Hallmark (barplot)
#### P1d: Hallmark activity and marker expression plots
##### Pancreas, Kidney, Colorectal, Breast, Brain
#### P1e: Moran's I bridge plots

### P2 
#### P2a2: ESTIMATE score in Colorectal tissue 
#### P2a3: ESTIMATE clusters in Colorectal tissue 
#### P2b: ESTIMATE validation boxplots
#### P2c: Pan-Cancer heatmap
#### P2d: hallmark expression example for Colorectal
#### P2e: Moran's I bridge plots for each compartment

### P3
#### P3a: Clonal maximum hallmark activity difference (heatmap)
#### P3b: Clonal clusters 
#### P3c: Hallmark activity across clones (boxplot)

### P4
#### P4b2: estimate clusters in kidney sample
#### P4b3: avoiding immune destruction in TME for kidney sample
#### P4b4: avoiding immune destruction radar in neoplastic for kidney sample
#### P4b5: avoiding immune destruction in TME and evading growth suppressors in neoplastic for kidney sample

### P5
#### P5a1: Feature importance to predict Cancer Hallmarks (boxplot)
#### P5a2: Feature importance to predict TME Hallmarks (boxplot)
#### P5b1: Stacked feature importance splited by spatial dependency in predicting Cancer Hallmarks (circos plot)
#### P5b2: Stacked feature importance splited by spatial dependency in predicting TME Hallmarks (circos plot)
#### P5c1: Random Forest for Cancer Hallmarks results (circos plot)
#### P5c2: Random Forest for Neoplastic Hallmarks results (circos plot)






# Code
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(tidyverse)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(ggforce)
library(scales)
#path to SamplesMetadata.R file
source("../utils/SamplesMetadata.R")
#path to PlottingMod.R file
source("../utils/PlottingMod.R")


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
H.genes.tbl <- dplyr::arrange(H.genes.tbl, n)
H.genes.tbl$Hallmark <- factor(H.genes.tbl$Hallmark, levels = H.genes.tbl$Hallmark)

ggbarplot(H.genes.tbl, x = "Hallmark", y = "n", fill = "black", color = "white", width = 0.9) +
  rotate_x_text(angle = -45, vjust = 2, hjust = 0) +  ylab("n.Genes") + coord_flip() + 
  theme(axis.ticks.y = element_blank())
ggsave("", bg="white")


#### P1d
#load spot_info ()
##### Pancreas
sample <- "Pancreas5"
#load enhanced object with imputed genes
sce <- readRDS("")
sce <- sce[,-1]
h <- "H11" #senescent cells
colData(sce)[,h] <- spot_info[spot_info$sample==sample, h]

v <- .make_triangle_subspots(colData(sce), fill = h)
p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c(option = "B") + labs(fill = "") + ggtitle(hallmark_names_list[[h]]) + theme(plot.title = element_text(size=35, hjust = 0.5),
                                                                                                   legend.position = c(1.05,0.3), legend.key.size = unit(1, "cm"))
#save file 
pdf("", width = 10, height = 10) 
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


##### Kidney
sample <- "Kidney2"
#load enhanced object with imputed genes
sce <- readRDS("")
sce <- sce[,-1]
h <- "H1" #Sustaining proliferative signaling
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
sample <- "Colorectal8"
#load enhanced object with imputed genes
sce <- readRDS("")
sce <- sce[,-1]
h <- "H8" #Genome Instability and Mutation
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
sample <- "Breast4"
#load enhanced object with imputed genes
sce <- readRDS("")
sce <- sce[,-1]
h <- "H12" #Nonmutational Epigenetic Reporgramming
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
sample <- "Glioblastoma1"
#load enhanced object with imputed genes
sce <- readRDS("")
sce <- sce[,-1]
h <- "H10" #Deregulating Cellular Energetics
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

####P1e
#load hallmark spatial autocorrelation for all samples
df <- read.table("", sep = "\t")

colnames(df)[2:14] <- hallmark_names
#ridge plot for each hallmark

df.long <- df %>% reshape2::melt(id.vars="sample")
df.long$variable <- factor(df.long$variable, levels = hallmark_names[rev(c(10,1,7,6,5,3,2,8,9,12,11,13,4))])
pdf("", height = 10, width = 8.5)
p <- ggplot(df.long, aes(y=variable, x=value, fill=variable)) + geom_density_ridges() + xlim(c(0,1)) +
  scale_fill_manual(values = do.call(c, color_codes))+ theme_ridges(center_axis_labels = T) + 
  geom_vline(xintercept = 0, linetype="dashed", col="gray32") + 
  theme(legend.position = "none", panel.grid.major.x = element_blank()) + labs(x="Moran's I", y="") 
plot(p)
dev.off()

### P2 
#image
sample <- "Colorectal1"
hires <- readRDS("")
#spot location
v <- readRDS("")
v <- v[v$spot != "subspot_1.1",]
#sub-spot data
sub_data <- read.table("")
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

#### P2b
#supplementary table S7
df <- readxl::read_xlsx("")
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

#### P2c
# run (PanCancerHeatmap.R)

#### P2d1
spot_info <- read.table("")
sample <- "Colorectal1"
print(sample)
hallmark <- "H13" #unlocking phenotypic plasticity
sce <- readRDS("")
sce <- sce[,-1]
colData(sce)[, hallmark] <- spot_info[spot_info$sample==sample, hallmark]

v <- .make_triangle_subspots(colData(sce), fill = hallmark)

ref_v <- readRDS("")
hires <- readRDS("")

for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}
pdf()
hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn(hallmark, colours = viridisLite::rocket(1000, alpha = 1, begin = 0, end = 1, direction = 1)[200:1000])
dev.off()
## violin
Colorectal1 <- spot_info %>% filter(sample == "Colorectal1") %>% 
  mutate(compartment = case_when(cluster %in% 1:2 ~ "Neoplastic",
                               cluster %in% 4:5 ~ "TME",
                               TRUE ~ "Buffer")) %>% 
  filter(compartment != "Buffer")

colors <- c("TME" = "lightgoldenrod1", "Neoplastic"= "orchid4")

ggplot(Colorectal1, aes(x = compartment, y = H13, color = compartment)) +
  geom_sina(size = 0.01, maxwidth = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, aes(color= NA)) +
  scale_color_manual(values = colors) +
  theme_classic() + ylab("Unlocking Phenotypic Plasticity")  +
  theme(axis.text=element_text(size=16, family = "Arial"),
        axis.title=element_text(size=18,face="bold")
  ) +
  stat_compare_means(family = "Arial", size = 5, label.y = 4) +
  geom_hline(yintercept=0, linetype= "dashed",
             color = "black", size=0.5)

#### P2d2
hallmark <- "H4" #enablign replicative immortality
sce <- readRDS("")
sce <- sce[,-1]
colData(sce)[, hallmark] <- spot_info[spot_info$sample==sample, hallmark]

v <- .make_triangle_subspots(colData(sce), fill = hallmark)

ref_v <- readRDS("")
hires <- readRDS("")

for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}
pdf()
hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn(hallmark, colours = viridisLite::rocket(1000, alpha = 1, begin = 0, end = 1, direction = 1)[200:1000])
dev.off()

## violin
ggplot(Colorectal1, aes(x = compartment, y = H4, color = compartment)) +
  geom_sina(size = 0.01, maxwidth = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, aes(color= NA)) +
  scale_color_manual(values = colors) +
  theme_classic() + ylab("Enabling Replicative Immortality")  +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold")
  ) +
  stat_compare_means(family = "Arial", size = 5, label.y = 4)+
  geom_hline(yintercept=0, linetype= "dashed",
             color = "black", size=0.5)


#### P2e
samples <- df$sample
df_all_samples <- c()
for (sample_id in samples) {
  df_whole <- df %>% filter(sample==sample_id) %>% remove_rownames %>% column_to_rownames(var="sample") %>% t()  %>% as.data.frame() 
  colnames(df_whole) <- c("cor")
  df_whole$tissue <- "Whole"
  df_whole$signature <- rownames(df_whole)
  
  
  #Neoplastic 
  df_neoplastic <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/moran/Neoplastic/", sample_id,"_HallmarkMoran.txt"), sep = "\t")
  df_neoplastic <- df_neoplastic %>% filter(!hallmark %in% paste0("H", c(1,3,5,6,7,9,13)))
  
  df_neoplastic <- df_neoplastic %>% select(-hallmark)
  df_neoplastic$tissue <- "Neoplastic"
  df_neoplastic$signature <- rownames(df_neoplastic)
  
  #TME 
  
  df_TME <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/moran/TME/", sample_id,"_HallmarkMoran.txt"), sep = "\t")
  df_TME <- df_TME %>% filter(hallmark %in% paste0("H", c(1,3,5,6,7,9,13)))
  
  df_TME <- df_TME %>% select(-hallmark)
  df_TME$tissue <- "TME"
  df_TME$signature <- rownames(df_TME)
  
  df_all <- rbind(df_whole, df_neoplastic, df_TME)
  df_all$id <- paste0(df_all$tissue, "_", df_all$signature2)
  df_all_samples <- rbind(df_all_samples, df_all)
}
#boxplot 
df_all_samples$tissue <- factor(df_all_samples$tissue, levels = c("Whole","Neoplastic", "TME"))
ggplot(df_all_samples, aes(x=tissue,y=cor, fill=tissue))  + 
  geom_violin() + geom_boxplot(width=0.1, outlier.size = 0.4, alpha=0.2, ) + ylim(c(0,1.05)) + 
  stat_compare_means(size=6, family="Times New Roman") + theme_classic() +
  labs(y="Moran's I", x="Tissue") + 
  ggtitle("Spatial autocorrelation of Hallmark activities") +
  scale_fill_manual(values = c("gray","orchid4", "lightgoldenrod1")) + 
  theme(axis.text = element_text(size=18), #axis.text.x = element_text( angle = -45, hjust = 0),
        axis.title = element_text(size=14, face="bold"),
        legend.position = "none", plot.title = element_text(size=18,hjust=0.5,face="bold"))
ggsave("", width = 7, height = 5)

### P3

#### P3a 
# load df_diff from CNVexperiment.R
col_fun = colorRamp2(c(0, max(df_diff)), c("white", "red"))
tissue_annot <- sapply(rownames(df_diff), function(x) {annotation_tissue[[x]]})
clonal_annot <- sapply(1:nrow(df_diff), function(x){ifelse(sum(df_diff[x, ]) == 0, "Monoclonal", "Polyclonal")})
col_tissue_map <- setNames(Seurat::DiscretePalette(15, palette = "polychrome") [c(1:3, 5, 6:11)], sort(unique(tissue_annot)))

d <- dist(df_diff[rowSums(df_diff) > 0,])
clus <- hclust(d, method = "complete")
plot(clus)
cutree(clus, 3)
diff <- c(cutree(clus, 3), setNames(rep(0,12),rownames(df_diff[rowSums(df_diff) == 0,])))
diff[diff==1] <- "Medium"
diff[diff==0] <- "Monoclonal"
diff[diff==2] <- "Low"
diff[diff==3] <- "High"
col_diff_map <- setNames(RColorBrewer::brewer.pal(12, "Paired")[c(1,3,5,7)], unique(diff))
row_ha = HeatmapAnnotation(`Hallmark activity difference` = factor(diff[rownames(df_diff)]),
                           `Tumor type` = tissue_annot,
                           col=list(`Tumor type` = col_tissue_map,
                                    `Hallmark activity difference` = col_diff_map), 
                           annotation_height = list(`Tumor type` = unit(2, "cm"),
                                                    `Clonality` = unit(1, "cm")))
palette <- do.call(c,color_codes[c(2,4,8,10,11,12)])
names(palette) <- NULL
top_ha = rowAnnotation(" " = anno_boxplot(df_diff, height = unit(2, "cm"),gp = gpar(fill = palette)))



pdf("", width = 16.5, height = 3)
Heatmap(t(df_diff), name = "Difference", border = "black", col=col_fun, show_row_names =T, 
        show_column_names =F,
        top_annotation =  row_ha, column_split = factor(clonal_annot, levels = c("Polyclonal", "Monoclonal")), cluster_column_slices=F,
        left_annotation = top_ha, 
        rect_gp = gpar(col = "gray", lwd = 0.1))
dev.off()

#### P3b 
cnv <- read.table("")
samples <- c("Breast7", "Colorectal8","Ovarian1", "Liver3")
cnv_clusters <- list(M3=0:2, Intestine=0:2, HCC2T=0:2, OV4A=0)
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
  hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol-500, group=~spot, fill=~as.factor(fill))) +  theme_void() +
    labs(fill="CNV clones") +  scale_fill_manual(values = c("#E69F00", "#0072B2", "#009E73"))
  dev.off()
}

#### P3c
samples <- c("Breast7", "Colorectal8","Ovarian1", "Liver3")
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

### P4
#### P4b2
#sample Kidney5
hires <- readRDS("")
#spot location
v <- readRDS("")
v <- v[v$spot != "subspot_1.1",]
#sub-spot data
sub_data <- read.table("")
cluster <- c("5" = "lightgoldenrod1", "4" = "lightgoldenrod3", "3" = "lightpink2", "2"= "orchid3", "1"= "orchid4")

v$fill <- sub_data[v$spot,"estimate"]
#save file 
pdf("", width = 10, height = 10) 
hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_gradientn(colours = rev(cluster), breaks=c(min(v$fill), max(v$fill)),labels=c("Cancer pure","TME pure")) + labs(fill="") + 
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.key.width = unit(118, "point"),
        legend.text = element_text(size = 20), plot.title = element_text(size = 45, hjust = 0.5)) +
  ggtitle("ESTIMATE score")
dev.off()
#### P4b3
spot_info <- read.table("")
sample <- "Kidney5"
h_tme <- "H3" #Avoiding Immune Destruction
sce <- readRDS("")
sce <- sce[,-1]
colData(sce)[, c(h_tme, "clusters")] <- spot_info[spot_info$sample==sample, c(h_tme, "clusters")]

v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2,3),], fill = "clusters")
v$fill <- "NA"
v2 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(4,5),], fill = h_tme)

ref_v <- readRDS("")
hires <- readRDS("")

for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}
for (spot in unique(v2$spot)) {
  v2$imagecol[v2$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v2$imagerow[v2$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}

pdf("", width = 7, height = 7)
hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_manual("Buffer", values = "gray41")  + 
  ggnewscale::new_scale("fill") + geom_polygon(data=v2,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_viridis_c()
dev.off()

#### P4b4
h_tme <- "H3_cancer" #Avoiding Immune Destruction (RADAR scores)
sce <- readRDS("")
sce <- sce[,-1]
#(TMERadar.R)
cancer_radar <- read.table("") 
colData(sce)[, c("clusters")] <- spot_info[spot_info$sample==sample, c("clusters")]
sce$H3_cancer <- NA
colData(sce)[sce$clusters %in% 1:2, h_tme] <- cancer_radar[cancer_radar$sample==sample, h_tme]

v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill = h_tme)
v2 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(3,4,5),], fill = "clusters")
v2$fill <- "NA"

ref_v <- readRDS("")
hires <- readRDS("")

for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}
for (spot in unique(v2$spot)) {
  v2$imagecol[v2$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v2$imagerow[v2$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}

pdf("", width = 7, height = 7)
hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_viridis_c()  + 
  ggnewscale::new_scale("fill") + geom_polygon(data=v2,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_manual("Buffer", values = "gray41")
dev.off()

#### P4b5
spot_info <- read.table("")
sample <- "Kidney5"
h_cancer <- "H2" #Evading Growth Suppressors
h_tme <- "H3" #Avoiding Immune Destruction
sce <- readRDS("")
sce <- sce[,-1]
colData(sce)[, c(h_cancer, h_tme, "clusters")] <- spot_info[spot_info$sample==sample, c(h_cancer, h_tme, "clusters")]

v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2,3),], fill = h_cancer)
v2 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(4,5),], fill = h_tme)

ref_v <- readRDS("")
hires <- readRDS("")

for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}

for (spot in unique(v2$spot)) {
  v2$imagecol[v2$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v2$imagerow[v2$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}

pdf("", width = 7, height = 7)
hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn(h_cancer, colours = viridisLite::rocket(1000, alpha = 1, begin = 0, end = 1, direction = 1)[200:1000]) + 
  ggnewscale::new_scale("fill") + geom_polygon(data=v2,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn(h_tme, colours = viridisLite::mako(1000, alpha = 1, begin = 0, end = 1, direction = 1)[200:1000]) 
dev.off()

### P5
#### P5a1
# load RF_results_Cancer from (CancerCircos.R)
# extract target, predictor and spatial direction and create long format df
direction_df <- melt(RF_results_Cancer[,c(14:20,13)], id.vars = "response")
colnames(direction_df)[3] <- "dir"
imp_df <- melt(RF_results_Cancer[,c(2:8,13)], id.vars = "response")
imp_df <- cbind(imp_df, direction_df[,3])
colnames(imp_df)[4] <- "dir"
imp_df$dir <- as.numeric(imp_df$dir)
imp_df$variable <- factor(imp_df$variable, levels = c("H3", "H1", "H9", "H6", "H7", "H13", "H5"))
imp_df <- arrange(imp_df, variable)

## combined with geom_sina
ggplot(imp_df, aes(x = variable, y = value, color = dir)) +
  geom_sina(size = 2, maxwidth = 0.8) +
  scale_color_gradient2(
    low = muted("blue"),
    mid = "white",
    high = muted("darkgreen"),
    midpoint = 0,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar"
  ) + geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  theme_classic() + ylab("Feature Importance Fraction")  +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold")
  ) +
  stat_compare_means(ref.group = c("H3"), size = 4.5, label = "p.format",family = "Arial", label.y = 0.9) +
  stat_compare_means(family = "Arial", size = 5, label.y = 1)

#### P5a2
# load RF_results_TME from (TMECircos.R)
# extract target, predictor and spatial direction and create long format df
direction_df <- melt(RF_results_TME[,c(13:18,12)], id.vars = "response")
colnames(direction_df)[3] <- "dir"
imp_df <- melt(RF_results_TME[,c(2:7,12)], id.vars = "response")
imp_df <- cbind(imp_df, direction_df[,3])
colnames(imp_df)[4] <- "dir"
imp_df$dir <- as.numeric(imp_df$dir)
imp_df$variable <- factor(imp_df$variable, levels = c("H10", "H11", "H2", "H4", "H12", "H8"))
imp_df <- arrange(imp_df, variable)

## combined with geom_sina
ggplot(imp_df, aes(x = variable, y = value, color = dir)) +
  geom_sina(size = 2, maxwidth = 0.8) +
  scale_color_gradient2(
    low = muted("blue"),
    mid = "white",
    high = muted("darkgreen"),
    midpoint = 0,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar"
  ) + geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  theme_classic() + ylab("Feature Importance Fraction")  +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold")
  ) +
  stat_compare_means(ref.group = c("H10"), size = 5, label = "p.format",family = "Arial", label.y = 0.9)+
  stat_compare_means(family = "Arial", size = 5, label.y = 1)

#### P5b1
RF_results_Cancer <- data.frame(RF_results_Cancer, row.names = 1)
RF_Cancer_Imp <- data.frame(t(RF_results_Cancer[1:7]))
RF_Cancer_dir <- data.frame(t(RF_results_Cancer[13:19]))
RF_Cancer_dir <- RF_Cancer_dir[match(paste0(rownames(RF_Cancer_Imp), "_cor"), rownames(RF_Cancer_dir)), ]
# sort each column and count how many hallmarks reach 0.75
RI_hallmarks_Cancer <- data.frame(matrix(nrow = 2436, ncol = 2)) # 348 * 7
colnames(RI_hallmarks_Cancer)[1:2] <- c("n_Hallmarks", "RI")
RI_hallmarks_Cancer$n_Hallmarks <- c(rep(1,348), rep(2,348), rep(3,348), rep(4,348), rep(5,348), rep(6,348), rep(7,348))
RI_hallmarks_Cancer$n_models <- c(rep(1:348,7))
RI_hallmarks_Cancer$sample <- ""
RI_hallmarks_Cancer$response <- ""
RI_hallmarks_Cancer$predictor <- ""
RI_hallmarks_Cancer$dir <- ""
for(h in 1:7){
  for (c in 1:ncol(RF_Cancer_Imp)){
    RF_Cancer_Imp <- arrange(RF_Cancer_Imp, desc(RF_Cancer_Imp[,c]))
    RF_Cancer_dir <- RF_Cancer_dir[match(paste0(rownames(RF_Cancer_Imp), "_cor"), rownames(RF_Cancer_dir)), ]
    RI_hallmarks_Cancer[RI_hallmarks_Cancer$n_models== c & RI_hallmarks_Cancer$n_Hallmarks == h,]$RI <- RF_Cancer_Imp[h,c]
    RI_hallmarks_Cancer[RI_hallmarks_Cancer$n_models== c & RI_hallmarks_Cancer$n_Hallmarks == h,]$dir <- RF_Cancer_dir[h,c]
    if (length(strsplit(names(RF_Cancer_Imp[c]), split = "_")[[1]]) < 3){
      RI_hallmarks_Cancer[RI_hallmarks_Cancer$n_models== c & RI_hallmarks_Cancer$n_Hallmarks == h,]$response <- strsplit(names(RF_Cancer_Imp[c]), "_")[[1]][2]
      RI_hallmarks_Cancer[RI_hallmarks_Cancer$n_models== c & RI_hallmarks_Cancer$n_Hallmarks == h,]$sample <- strsplit(names(RF_Cancer_Imp[c]), "_")[[1]][1]
    }
    else {
      RI_hallmarks_Cancer[RI_hallmarks_Cancer$n_models== c & RI_hallmarks_Cancer$n_Hallmarks == h,]$response <- strsplit(names(RF_Cancer_Imp[c]), "_")[[1]][3]
      RI_hallmarks_Cancer[RI_hallmarks_Cancer$n_models== c & RI_hallmarks_Cancer$n_Hallmarks == h,]$sample <- paste0(strsplit(names(RF_Cancer_Imp[c]), split = "_")[[1]][1],"_",strsplit(names(RF_Cancer_Imp[c]), split = "_")[[1]][2])
    }
    RI_hallmarks_Cancer[RI_hallmarks_Cancer$n_models== c & RI_hallmarks_Cancer$n_Hallmarks == h,]$predictor <- rownames(RF_Cancer_Imp[h,])
  }
}
type <- RF_results_Cancer[c("sample", "type")]
type <- type[duplicated(type) == F,]
RI_hallmarks_Cancer$type <- ""
for (s in type$sample) {
  RI_hallmarks_Cancer[RI_hallmarks_Cancer$sample == s,]$type <- type[type$sample == s,]$type
}
RI_hallmarks_Cancer$dir <- as.numeric(RI_hallmarks_Cancer$dir)
RI_hallmarks_Cancer$dir_cat <- ""
for (r in 1:nrow(RI_hallmarks_Cancer)) {
  if (RI_hallmarks_Cancer[r,]$dir >= 0.6){
    RI_hallmarks_Cancer[r,]$dir_cat <- "Positive"
  }
  else if (RI_hallmarks_Cancer[r,]$dir <= -0.6){
    RI_hallmarks_Cancer[r,]$dir_cat <- "Negative"
  } else {RI_hallmarks_Cancer[r,]$dir_cat <- "nonlinear"}
}
## fraction for rank 1
color_codes <-  c("H1"= "#15CE59",
                  "H2" = "#701717",
                  "H3" = "#CB3BBD",
                  "H4" = "#4969D4",
                  "H5" = "#E6880D",
                  "H6" = "#000000",
                  "H7" = "#EE0F16",
                  "H8" = "#132892",
                  "H9" = "#8E909B",
                  "H10" = "#71189E",
                  "H11" = "#05F3EB",
                  "H12" = "#890269",
                  "H13" = "#95641A")
RI_hallmarks_Cancer_fr <- RI_hallmarks_Cancer

RI_hallmarks_Cancer_fr <- filter(RI_hallmarks_Cancer, n_Hallmarks == 1)
# barplot
RI_hallmarks_Cancer_fr$dir_cat <- factor(RI_hallmarks_Cancer_fr$dir_cat, levels = c("Positive", "Negative", "nonlinear"))
RI_hallmarks_Cancer_fr <- arrange(RI_hallmarks_Cancer_fr, dir_cat)
RI_hallmarks_Cancer_fr$predictor <- factor(RI_hallmarks_Cancer_fr$predictor, levels = c("H3", "H1", "H9", "H6", "H7", "H13", "H5"))


ggplot (RI_hallmarks_Cancer_fr,aes (x = response, fill = predictor, y = RI)) + geom_bar(stat = "identity", position = "stack") + 
  facet_wrap (~dir_cat) + scale_fill_manual(values = color_codes[c("H3", "H1", "H9", "H6", "H7", "H13", "H5")]) +
  theme_classic()
#### P5b2
RF_results_TME <- data.frame(RF_results_TME, row.names = 1)

RF_TME_Imp <- data.frame(t(RF_results_TME[1:6]))
View(RF_TME_Imp)

RF_TME_dir <- data.frame(t(RF_results_TME[12:17]))
RF_TME_dir <- RF_TME_dir[match(paste0(rownames(RF_TME_Imp), "_cor"), rownames(RF_TME_dir)), ]
View(RF_TME_dir)

#
RI_hallmarks_TME <- data.frame(matrix(nrow = 2436, ncol = 2))
colnames(RI_hallmarks_TME)[1:2] <- c("n_Hallmarks", "RI")
RI_hallmarks_TME$n_Hallmarks <- c(rep(1,406), rep(2,406), rep(3,406), rep(4,406), rep(5,406), rep(6,406))
RI_hallmarks_TME$n_models <- c(rep(1:406,6))
RI_hallmarks_TME$sample <- ""
RI_hallmarks_TME$response <- ""
RI_hallmarks_TME$predictor <- ""
RI_hallmarks_TME$dir <- ""

for(h in 1:6){
  for (c in 1:ncol(RF_TME_Imp)){
    RF_TME_Imp <- arrange(RF_TME_Imp, desc(RF_TME_Imp[,c]))
    RF_TME_dir <- RF_TME_dir[match(paste0(rownames(RF_TME_Imp), "_cor"), rownames(RF_TME_dir)), ]
    RI_hallmarks_TME[RI_hallmarks_TME$n_models== c & RI_hallmarks_TME$n_Hallmarks == h,]$RI <- RF_TME_Imp[h,c]
    RI_hallmarks_TME[RI_hallmarks_TME$n_models== c & RI_hallmarks_TME$n_Hallmarks == h,]$dir <- RF_TME_dir[h,c]
    
    if (length(strsplit(names(RF_TME_Imp[c]), split = "_")[[1]]) < 3){
      RI_hallmarks_TME[RI_hallmarks_TME$n_models== c & RI_hallmarks_TME$n_Hallmarks == h,]$response <- strsplit(names(RF_TME_Imp[c]), "_")[[1]][2]
      RI_hallmarks_TME[RI_hallmarks_TME$n_models== c & RI_hallmarks_TME$n_Hallmarks == h,]$sample <- strsplit(names(RF_TME_Imp[c]), "_")[[1]][1]
    }
    else {
      RI_hallmarks_TME[RI_hallmarks_TME$n_models== c & RI_hallmarks_TME$n_Hallmarks == h,]$response <- strsplit(names(RF_TME_Imp[c]), "_")[[1]][3]
      RI_hallmarks_TME[RI_hallmarks_TME$n_models== c & RI_hallmarks_TME$n_Hallmarks == h,]$sample <- paste0(strsplit(names(RF_TME_Imp[c]), split = "_")[[1]][1],"_",strsplit(names(RF_TME_Imp[c]), split = "_")[[1]][2])
    }
    RI_hallmarks_TME[RI_hallmarks_TME$n_models== c & RI_hallmarks_TME$n_Hallmarks == h,]$predictor <- rownames(RF_TME_Imp[h,])
  }
}

type <- RF_results_TME[c("sample", "type")]
type <- type[duplicated(type) == F,]

RI_hallmarks_TME$type <- ""
for (s in type$sample) {
  RI_hallmarks_TME[RI_hallmarks_TME$sample == s,]$type <- type[type$sample == s,]$type
}

RI_hallmarks_TME$dir <- as.numeric(RI_hallmarks_TME$dir)
RI_hallmarks_TME$dir_cat <- ""
for (r in 1:nrow(RI_hallmarks_TME)) {
  if (RI_hallmarks_TME[r,]$dir >= 0.6){
    RI_hallmarks_TME[r,]$dir_cat <- "Positive"
  }
  else if (RI_hallmarks_TME[r,]$dir <= -0.6){
    RI_hallmarks_TME[r,]$dir_cat <- "Negative"
  } else {RI_hallmarks_TME[r,]$dir_cat <- "nonlinear"}
  
}


RI_hallmarks_TME_fr <- RI_hallmarks_TME

RI_hallmarks_TME_fr <- filter(RI_hallmarks_TME, n_Hallmarks == 1)
# barplot
RI_hallmarks_TME_fr$dir_cat <- factor(RI_hallmarks_TME_fr$dir_cat, levels = c("Positive", "Negative", "nonlinear"))
RI_hallmarks_TME_fr <- arrange(RI_hallmarks_TME_fr, dir_cat)
RI_hallmarks_TME_fr$predictor <- factor(RI_hallmarks_TME_fr$predictor, levels = c("H10", "H11", "H2", "H4", "H12", "H8"))


ggplot (RI_hallmarks_TME_fr,aes (x = response, fill = predictor, y = RI)) + geom_bar(stat = "identity", position = "stack") + 
  facet_wrap (~dir_cat) + scale_fill_manual(values = color_codes[c("H10", "H11", "H2", "H4", "H12", "H8")]) +
  theme_classic()

#### P5c1
# (CancerCircos.R)
#### P5c2
# (TMECircos.R)



