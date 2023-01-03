################################################################################
#                      Sample: 4A
################################################################################

library(STutility)

infoTable <- data.frame(samples="Desktop/IJC/datasets/IGTP/4A/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5")
infoTable$spotfiles="Desktop/IJC/datasets/IGTP/4A/spatial/tissue_positions_list.csv"
infoTable$imgs="Desktop/IJC/datasets/IGTP/4A/spatial/tissue_hires_image.png"
infoTable$json="Desktop/IJC/datasets/IGTP/4A/spatial/scalefactors_json.json"

OV4A <- InputFromTable(infotable = infoTable, 
                       platform =  "Visium",
                       minUMICountsPerGene = 5)


protein_genes <- annotLookup$gene_name[annotLookup$gene_type %in% c("protein_coding", "TR_V_gene", "TR_D_gene", "TR_J_gene", "TR_C_gene", "IG_LV_gene", "IG_V_gene", "IG_J_gene", "IG_C_gene" , "IG_D_gene")]
selected_genes <- rownames(OV4A)
selected_genes <- protein_genes[protein_genes %in% selected_genes]
selected_genes <- selected_genes[!grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", selected_genes)]
selected_genes <- selected_genes[!grepl("^MT-", selected_genes)]

OV4A <- SubsetSTData(OV4A, features = selected_genes)

# Load images into a Seurat objectm, mask and plot images
OV4A <- LoadImages(OV4A, verbose = TRUE) %>% MaskImages()
ImagePlot(OV4A)

img <- Seurat::Read10X_Image(image.dir = 'Desktop/IJC/datasets/IGTP/4A/images/')
Seurat::DefaultAssay(object = img) <- 'RNA'
rownames(img@coordinates) <- paste0(rownames(img@coordinates), "_1")
img <- img[colnames(x = OV4A)]
OV4A[['OC']] <- img


OV4A <- ManualAnnotation(OV4A)
OV4A <- SubsetSTData(OV4A, labels == "Default")

OV4A <- SCTransform(OV4A, return.only.var.genes = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.1, assay = "RNA")
VlnPlot(OV4A, features = c("nCount_RNA", "nFeature_RNA"))




OV4A <- RunNMF(OV4A, assay = "SCT", nfactors = 15)
OV4A <- AddMetaData(OV4A, as.data.frame(OV4A@reductions$NMF@cell.embeddings))

p <- SpatialFeaturePlot(OV4A, features = colnames(OV4A@meta.data)[8:22], ncol=3)
pdf("Desktop/IJC/datasets/IGTP/4A/figures/factors15/factors15.pdf", width = 10, height = 30)
plot(p)
dev.off()

for (i in 1:15) {
  file_name = paste0("Desktop/IJC/datasets/IGTP/4A/figures/factors15/loadings/loadings_factor", i, ".jpeg")
  ggsave(file_name, plot = FactorGeneLoadingPlot(OV4A, factor = i, topn = 30), bg = "white")
}



################################################################################
#                      Sample: 5A
################################################################################



infoTable <- data.frame(samples="Desktop/IJC/datasets/IGTP/5A/filtered_feature_bc_matrix.h5")
infoTable$spotfiles="Desktop/IJC/datasets/IGTP/5A/images/tissue_positions_list.csv"
infoTable$imgs="Desktop/IJC/datasets/IGTP/5A/images/tissue_hires_image.png"
infoTable$json="Desktop/IJC/datasets/IGTP/5A/images/scalefactors_json.json"

OV5A <- InputFromTable(infotable = infoTable, 
                       platform =  "Visium",
                       minUMICountsPerGene = 5)


protein_genes <- annotLookup$gene_name[annotLookup$gene_type %in% c("protein_coding", "TR_V_gene", "TR_D_gene", "TR_J_gene", "TR_C_gene", "IG_LV_gene", "IG_V_gene", "IG_J_gene", "IG_C_gene" , "IG_D_gene")]
selected_genes <- rownames(OV5A)
selected_genes <- protein_genes[protein_genes %in% selected_genes]
selected_genes <- selected_genes[!grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", selected_genes)]
selected_genes <- selected_genes[!grepl("^MT-", selected_genes)]

OV5A <- SubsetSTData(OV5A, features = selected_genes)

# Load images into a Seurat objectm, mask and plot images
OV5A <- LoadImages(OV5A, verbose = TRUE) %>% MaskImages()
ImagePlot(OV5A)

img <- Seurat::Read10X_Image(image.dir = 'Desktop/IJC/datasets/IGTP/5A/images/')
Seurat::DefaultAssay(object = img) <- 'RNA'
rownames(img@coordinates) <- paste0(rownames(img@coordinates), "_1")
img <- img[colnames(x = OV5A)]
OV5A[['OC']] <- img


OV5A <- ManualAnnotation(OV5A)
OV5A <- SubsetSTData(OV5A, labels == "Default")

OV5A <- SCTransform(OV5A, return.only.var.genes = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.1, assay = "RNA")
VlnPlot(OV5A, features = c("nCount_RNA", "nFeature_RNA"))
OV5A <- subset(OV5A, nFeature_RNA > 1500)






OV5A <- RunNMF(OV5A, assay = "SCT", nfactors = 15)
OV5A <- AddMetaData(OV5A, as.data.frame(OV5A@reductions$NMF@cell.embeddings))

p <- SpatialFeaturePlot(OV5A, features = colnames(OV5A@meta.data)[8:22], ncol=3)
pdf("Desktop/IJC/datasets/IGTP/5A/figures/factors15/factors15.pdf", width = 10, height = 30)
plot(p)
dev.off()

for (i in 1:15) {
  file_name = paste0("Desktop/IJC/datasets/IGTP/5A/figures/factors15/loadings/loadings_factor", i, ".jpeg")
  ggsave(file_name, plot = FactorGeneLoadingPlot(OV5A, factor = i, topn = 30), bg = "white")
}

################################################################################
#                      Sample: 5B
################################################################################



infoTable <- data.frame(samples="Desktop/IJC/datasets/IGTP/5B/filtered_feature_bc_matrix.h5")
infoTable$spotfiles="Desktop/IJC/datasets/IGTP/5B/images/tissue_positions_list.csv"
infoTable$imgs="Desktop/IJC/datasets/IGTP/5B/images/tissue_hires_image.png"
infoTable$json="Desktop/IJC/datasets/IGTP/5B/images/scalefactors_json.json"

OV5B <- InputFromTable(infotable = infoTable, 
                       platform =  "Visium",
                       minUMICountsPerGene = 5)


protein_genes <- annotLookup$gene_name[annotLookup$gene_type %in% c("protein_coding", "TR_V_gene", "TR_D_gene", "TR_J_gene", "TR_C_gene", "IG_LV_gene", "IG_V_gene", "IG_J_gene", "IG_C_gene" , "IG_D_gene")]
selected_genes <- rownames(OV5B)
selected_genes <- protein_genes[protein_genes %in% selected_genes]
selected_genes <- selected_genes[!grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", selected_genes)]
selected_genes <- selected_genes[!grepl("^MT-", selected_genes)]

OV5B <- SubsetSTData(OV5B, features = selected_genes)

# Load images into a Seurat objectm, mask and plot images
OV5B <- LoadImages(OV5B, verbose = TRUE) %>% MaskImages()
ImagePlot(OV5B)

img <- Seurat::Read10X_Image(image.dir = 'Desktop/IJC/datasets/IGTP/5B/images/')
Seurat::DefaultAssay(object = img) <- 'RNA'
rownames(img@coordinates) <- paste0(rownames(img@coordinates), "_1")
img <- img[colnames(x = OV5B)]
OV5B[['OC']] <- img


OV5B <- ManualAnnotation(OV5B)
OV5B <- SubsetSTData(OV5B, labels == "Default")

VlnPlot(OV5B, features = c("nCount_RNA", "nFeature_RNA"))
OV5B <- subset(OV5B, nFeature_RNA > 1500)
OV5B <- SCTransform(OV5B, return.only.var.genes = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.1, assay = "RNA")



OV5B <- RunNMF(OV5B, assay = "SCT", nfactors = 15)
OV5B <- AddMetaData(OV5B, as.data.frame(OV5B@reductions$NMF@cell.embeddings))

p <- SpatialFeaturePlot(OV5B, features = colnames(OV5B@meta.data)[8:22], ncol=3)
pdf("Desktop/IJC/datasets/IGTP/5B/figures/factors15/factors15.pdf", width = 10, height = 30)
plot(p)
dev.off()

for (i in 1:15) {
  file_name = paste0("Desktop/IJC/datasets/IGTP/5B/figures/factors15/loadings/loadings_factor", i, ".jpeg")
  ggsave(file_name, plot = FactorGeneLoadingPlot(OV5B, factor = i, topn = 30), bg = "white")
}

