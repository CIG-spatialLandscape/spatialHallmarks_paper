OV4A$patient <- "4"
OV5A$patient <- "5"
OV5B$patient <- "5"
OV <- merge(x = OV4A, y = list(OV5A, OV5B))

OV <- SCTransform(OV,return.only.var.genes = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.1, assay = "RNA", vars.to.regress = "patient")

OV <- OV4A


library(pheatmap)


names(sort(OV@reductions$NMF@feature.loadings[, 11], decreasing = T)[1:10])


genes <- c()
for (factor in 1:15) {
  genes <- c(genes, names(sort(OV4A@reductions$NMF@feature.loadings[, factor], decreasing = T)[1:5]))
}

mat <- t(scale(OV4A@reductions$NMF@feature.loadings[genes,]))

pheatmap(mat, 
         cluster_rows = F, cluster_cols = F, 
         color = rev(viridis::viridis(100,option = "magma")),
          legend_labels  = c("low", "high"), legend_breaks = c(min(mat),max(mat)),border_color = NA
        )


palette <- RColorBrewer::brewer.pal(10, name = "Paired")

p <- SpatialDimPlot(OV, group.by = "bayes.cluster", cols = palette) + 
  theme(legend.key = element_blank(), legend.direction = "horizontal", legend.position = "top", legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 4) ))
ggsave("Desktop/IJC/datasets/IGTP/figuresPaper/4A/BayesSpace/clusters.png", width = 7, height = 7)

p <- SpatialDimPlot(OV5A, group.by = "bayes.cluster2", cols = palette) + 
  theme(legend.key = element_blank(), legend.direction = "horizontal", legend.position = "top", legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 4) ))
ggsave("Desktop/IJC/datasets/IGTP/figuresPaper/4A/BayesSpace/clusters.png", width = 7, height = 7)

p <- SpatialDimPlot(OV5B, group.by = "bayes.cluster", cols = palette) + 
  theme(legend.key = element_blank(), legend.direction = "horizontal", legend.position = "top", legend.title = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(size = 4) ))
ggsave("Desktop/IJC/datasets/IGTP/figuresPaper/4A/BayesSpace/clusters.png", width = 7, height = 7)