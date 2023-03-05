##################################################
## Project: Hallmarks paper
## Script purpose: Supplementary Figures
## Date: 23/1/2023
## Author: Sergi Cervilla* & Mustafa Sibai*
##################################################

library(tidyverse)
library(dplyr)
library(ggplot2)
library(scales)
source("Desktop/IJC/git/spatialHallmarks_paper/utils/SamplesMetadata.R")
# S1 - number of samples per cancer type / Embedding method
df <- read.table("Downloads/RF_Imp_direction_cancer.txt")
df <-  df %>% dplyr::select(type, sample) %>% unique()
df$method <- sapply(df$sample, function(x) {annotation_method[[x]]})
df_plot <- df

df_plot <- df %>% group_by(type) %>% dplyr::count()

ggplot(df, aes(x=type, fill=method)) + geom_bar() + theme_classic() + 
  scale_y_continuous(breaks= pretty_breaks()) + 
  labs(x="Tumor type", y="Number of samples", fill="Embedding Method") + 
  scale_fill_manual(values = Seurat::DiscretePalette(15, palette = "polychrome")[c(11, 14)])

# S2 - number of samples per cancer type / sample origin method
df <- read.table("Downloads/RF_Imp_direction_cancer.txt")
df <-  df %>% dplyr::select(type, sample) %>% unique()
df$generated <- sapply(df$sample, function(x) {annotation_generated[[x]]})
df_plot <- df

ggplot(df_plot, aes(x=type, fill=generated)) + geom_bar() + theme_classic() + 
  scale_y_continuous(breaks= pretty_breaks()) + 
  labs(x="Tumor type", y="Number of samples", fill="Sample Origin") + 
  scale_fill_manual(values = c("darkorchid2", "salmon1"))

#S3 hallmarks per cluster per cancer type
df_hallmark <- read.table("Desktop/df_cnv.txt")
df1 <- df_hallmark %>% dplyr::select(sample, clusters, H1:H13) %>% group_by(sample, clusters) %>% summarise_all(mean)
df1$clusters <- factor(df1$clusters)

for (h in paste0("H", 1:13)) {
  print(h)
  ggplot(df1, aes_string(x="clusters", y=h, fill="clusters")) + geom_violin(alpha=0.6) + geom_boxplot(width=0.2) +
    geom_jitter(width = 0.2, size=0.3) + theme_classic() + 
    scale_fill_manual(values = c("5" = "lightgoldenrod1", "4" = "lightgoldenrod3", "3" = "lightpink2", "2"= "orchid3", "1"= "orchid4"))+
    labs(x="ESTIMATE clusters", y="Average Hallmark activity") + 
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 0, linetype ="dashed") + ggtitle(hallmark_names_list[[h]])
  ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/SuppFigures/", h, "_avgBoxplot.png"), width = 7, height = 7, bg = "white")
}

df_hallmark <- read.table("Desktop/df_cnv.txt")
df1 <- df_hallmark %>% dplyr::select(sample, clusters, H1:H13) %>% group_by(sample, clusters) %>% summarise_all(mean)
df1$clusters <- factor(df1$clusters)
df1$type <- sapply(df1$sample, function(x) {annotation_tissue[[x]]})
for (h in paste0("H", 1:13)) {
  print(h)
  ggplot(df1, aes_string(x="clusters", y=h, fill="clusters")) + geom_violin(alpha=0.6) + geom_boxplot(width=0.2) +
    geom_jitter(width = 0.2, size=0.3) + theme_classic() + 
    scale_fill_manual(values = c("5" = "lightgoldenrod1", "4" = "lightgoldenrod3", "3" = "lightpink2", "2"= "orchid3", "1"= "orchid4"))+
    labs(x="ESTIMATE clusters", y="Average Hallmark activity") + 
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 25),
          axis.text = element_text(size = 15), axis.title = element_text(size=20),strip.text = element_text(size=15)
          ) +
    geom_hline(yintercept = 0, linetype ="dashed") + ggtitle(hallmark_names_list[[h]]) + facet_wrap(~type, nrow = 1)
  ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/SuppFigures/", h, "_avgBoxplotCTP.png"), width = 20, height = 4, bg = "white")
}

## Correlation Plots
files <- list.files("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/output_df", full.names = F)
files <- stringr::str_remove(files, pattern = ".txt")


correlations <- list()
for (sample in files) {
  hallmarks <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/objects_mts/", sample, "_hallmarks.txt"), sep = "\t", header = T)
  hallmarks <- hallmarks[rownames(hallmarks)!="subspot_1.1",]
  correlations[[sample]] <- cor(hallmarks)
}

avg <- as.data.frame(do.call(rbind, correlations))
avg$h_names <- rep(paste0("H", 1:13), 58)
avg <- avg %>% group_by(h_names) %>% summarise(H1 = mean(H1),
                                               H2 = mean(H2),
                                               H3 = mean(H3),
                                               H4 = mean(H4),
                                               H5 = mean(H5),
                                               H6 = mean(H6),
                                               H7 = mean(H7),
                                               H8 = mean(H8),
                                               H9 = mean(H9),
                                               H10 = mean(H10),
                                               H11 = mean(H11),
                                               H12 = mean(H12),
                                               H13 = mean(H13)) %>% as.data.frame()
rownames(avg) <- avg$h_names
avg$h_names <- NULL

source("Desktop/IJC/git/spatialHallmarks_paper/utils/SamplesMetadata.R")
colnames(avg) <- sapply(colnames(avg), function(x){hallmark_names_list[[x]]})
rownames(avg) <- sapply(rownames(avg), function(x){hallmark_names_list[[x]]})

library(pheatmap)
library(paletteer)
png(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/correlation_plots/Pancancer_cor.png"), width = 350, height = 350, )
pheatmap(avg, scale = "none", color = paletteer_c("grDevices::Blue-Red 3", 41), 
         breaks = seq(-1, 1, by = 0.05), angle_col = 315, main = "Pan-Cancer", fontsize = 12
)
dev.off()

## S: Hallmark correlations

tumor_type <- list()
tumor_type[["Prostate"]] <- c("Acinar", "IC", "PC1", "PC2")
tumor_type[["Breast"]] <- c("Breast", "BreastA", "Ductal", "DuctalFFPE", "TNBCA", "M1", "M2", "M3", "M4")
tumor_type[["Kidney"]]  <- c("C20", "C21", "C34", "C51", "C7")
tumor_type[["Liver"]] <- c("cHC1T", "HCC1T", "HCC2T", "HCC5D", "ICC1L")
tumor_type[["Colorectal"]] <- c("Colorectal", "CRC1", "CRC2", "Intestine", "Co1", "Co2", "Co3", "Co4")
tumor_type[["Bladder"]] <- c("DU2", "DU3", "DU12", "DU8", "DU13")
tumor_type[["Ovarian"]] <- c("OV4A", "OVD1", "Ovarian", "CUP295", "OVFFPE")
tumor_type[["Pancreas"]] <- c("P259_H2A2", "P264", "P270", "P288", "P306")
tumor_type[["Glioblastoma"]] <- c("Glioblastoma", "UKF242T", "UKF260T", "UKF269T", "UKF275T")
tumor_type[["Lung"]] <- c("P1","P3", "P4", "P5", "P6", "P7","P8")


# Cancer-type
for (tumor in names(tumor_type)) {
  correlations <- list()
  for (sample in tumor_type[[tumor]]) {
    hallmarks <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/objects_mts/", sample, "_hallmarks.txt"), sep = "\t", header = T)
    hallmarks <- hallmarks[rownames(hallmarks)!="subspot_1.1",]
    correlations[[sample]] <- cor(hallmarks)
  }
  avg <- as.data.frame(do.call(rbind, correlations))
  avg$h_names <- rep(paste0("H", 1:13), length(tumor_type[[tumor]]))
  avg <- avg %>% group_by(h_names) %>% summarise(H1 = mean(H1),
                                                 H2 = mean(H2),
                                                 H3 = mean(H3),
                                                 H4 = mean(H4),
                                                 H5 = mean(H5),
                                                 H6 = mean(H6),
                                                 H7 = mean(H7),
                                                 H8 = mean(H8),
                                                 H9 = mean(H9),
                                                 H10 = mean(H10),
                                                 H11 = mean(H11),
                                                 H12 = mean(H12),
                                                 H13 = mean(H13)) %>% as.data.frame()
  rownames(avg) <- avg$h_names
  avg$h_names <- NULL
  png(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/correlation_plots/", tumor, "_cor.png"), width = 350, height = 350, )
  pheatmap(avg, scale = "none", color = paletteer_c("grDevices::Blue-Red 3", 41), 
           breaks = seq(-1, 1, by = 0.05), angle_col = 315, main = tumor, fontsize = 12
  )
  dev.off()
}

# SCD across cancer types
df <- read.table("Desktop/df_meta.txt")
df
df$type <- sapply(rownames(df), function(x){annotation_tissue[[x]]})
ggplot(df, aes(x=type, col=type, y=SCD)) + geom_boxplot() + geom_point() + 
  scale_color_manual(values = Seurat::DiscretePalette(15, palette = "polychrome") [c(1:3, 5, 6:10, 12)]) +
  theme_classic()


## S: 
