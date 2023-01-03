library(Seurat)
library(maple)
library(patchwork)
library(ggplot2)



OV4 <- readRDS("Desktop/IJC/datasets/IGTP/4A/RDS/OV4A.rds")
OV4$patient <- "4"
substr(colnames(OV4), start = nchar(colnames(OV4)[1]), stop =nchar(colnames(OV4)[1]))

OV5A <- readRDS("Desktop/IJC/datasets/IGTP/5A/RDS/OV5A.rds")
OV5A$patient <- "5"
OV5B <- readRDS("Desktop/IJC/datasets/IGTP/5B/RDS/OV5B.rds")
OV5B$patient <- "5"
colnames(OV5B) <- paste0(substr(colnames(OV5B), start = 0, stop =nchar(colnames(OV5B)[1])-1), "3") 
all <- merge(OV4, OV5A)
all <- merge(all, OV5B)


DefaultAssay(all) <- "SCT"
VariableFeatures(all) <- c(VariableFeatures(OV4),
                             VariableFeatures(OV5A),
                             VariableFeatures(OV5B))
all <- RunPCA(all)
fit_PCs_K10 <- fit_maple(all,K = 5,emb = "harmony",covars = "patient",nsim = 1000, burn = 500)


all$maple_labels_PCs_K10 <- factor(fit_PCs_K10$z,
                                     levels = sort(as.numeric(unique(fit_PCs_K10$z_init))))
Idents(all) <- "maple_labels_PCs_K10"
SpatialDimPlot(all,
                                   label = TRUE,
                                   combine = T,
                                   label.size = 4,
                                   repel = TRUE)
  z_plots_K10 <- z_plots_K10_list[[1]] + z_plots_K10_list[[3]] + 
  z_plots_K10_list[[2]] + z_plots_K10_list[[4]] + 
  plot_layout(nrow = 2, byrow = TRUE)
z_plots_K10



delta_plots_K10 <- spruce::plot_deltas(fit_PCs_K10)
delta_plots_K10


post_scores <- spruce::get_scores(fit_PCs_K10)
all@meta.data <- cbind(all@meta.data,post_scores)

# plot scores
SpatialFeaturePlot(all,features = "u_score",combine = T)
u_plots_K10 <- u_plots_K10_list[[1]] + u_plots_K10_list[[3]] + 
  u_plots_K10_list[[2]] + u_plots_K10_list[[4]] + 
  plot_layout(nrow = 2, byrow = TRUE)
u_plots_K10





k1_plots_K10_list <- SpatialFeaturePlot(all,features = "k1_score", combine = FALSE)
k1_plots_K10 <- k1_plots_K10_list[[1]] + k1_plots_K10_list[[3]] + 
  k1_plots_K10_list[[2]] + k1_plots_K10_list[[4]] + 
  plot_layout(nrow = 2, byrow = TRUE)
k1_plots_K10


# plot proportions
colnames(fit_PCs_K10$W) <- c("4","5")
fit_PCs_K10$z <- factor(fit_PCs_K10$z,
                              levels = sort(as.numeric(unique(fit_PCs_K10$z))))
fit_PCs_K10$W[,"5"] <- ifelse(fit_PCs_K10$W[,"5"] == 1, "5","4")

grouped_bar_plot_K10 <- plot_props_grouped(fit_PCs_K10, group = "5")+ 
  theme(text = element_text(size = 20))
grouped_bar_plot_K10



# plot alluvial
plot_props_alluvial(fit_PCs_K10, group = "Posterior")

pie_plot_K10 <- plot_props_pie(fit_PCs_K10, group = "Posterior") + 
  theme(text = element_text(size = 20))
pie_plot_K10



all <- RunHarmony(all, "patient", verbose = F)
all <- RunUMAP(all, reduction = "harmony", dims = 1:50)
DimPlot(all, reduction = "umap", group.by = "patient")




OV <- cbind(OV_4A, OV_5A, OV_5B, deparse.level = 1)


merge(OV_4A, OV_5B)

cbind(OV_5A, OV_5B)

l <- list("1"=OV_4A, "2"=OV_5A, "3"=OV_5B)
all_sce <- sce_cbind(l, exprs = "SCT", method = "union")
