OV4A <- RunNMF(OV4A, assay = "SCT", nfactors = 15, n.cores = 10)
OV4A <- AddMetaData(OV4A, metadata = as.data.frame(OV4A@reductions$NMF@cell.embeddings))

p <- SpatialFeaturePlot(OV4A, features = colnames(OV4A@meta.data)[10:24], ncol=3)
pdf("Desktop/IJC/datasets/IGTP/4A/figures/factors15_filtered/factors15.pdf", width = 10, height = 30)
plot(p)
dev.off()

for (i in 1:15) {
  file_name = paste0("Desktop/IJC/datasets/IGTP/4A/figures/factors15_filtered/loadings/loadings_factor", i, ".jpeg")
  ggsave(file_name, plot = FactorGeneLoadingPlot(OV4A, factor = i, topn = 30))
}


OV <- RunNMF(OV, assay = "SCT", nfactors = 15, n.cores = 10)
OV <- AddMetaData(OV, metadata = as.data.frame(OV@reductions$NMF@cell.embeddings))

p <- SpatialFeaturePlot(OV, features = colnames(OV@meta.data)[10:24], ncol=3)
pdf("Desktop/IJC/datasets/IGTP/all/figures/factors15_filtered/factors15.pdf", width = 10, height = 80)
plot(p)
dev.off()

for (i in 1:15) {
  file_name = paste0("Desktop/IJC/datasets/IGTP/all/figures/factors15_filtered/loadings/loadings_factor", i, ".jpeg")
  ggsave(file_name, plot = FactorGeneLoadingPlot(OV, factor = i, topn = 30))
}
