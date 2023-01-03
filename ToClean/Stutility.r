
infoTable <- data.frame(samples="Desktop/IJC/datasets/IGTP/4A/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5")
infoTable$spotfiles="Desktop/IJC/datasets/IGTP/4A/spatial/tissue_positions_list.csv"
infoTable$imgs="Desktop/IJC/datasets/IGTP/4A/spatial/tissue_hires_image.png"
infoTable$json="Desktop/IJC/datasets/IGTP/4A/spatial/scalefactors_json.json"



gene_list <- read.table("Desktop/annotLookup.xls", header = T)

protein_genes <- gene_list$external_gene_name[gene_list$gene_biotype %in% c("protein_coding", "TR_V_gene", "TR_D_gene", "TR_J_gene", "TR_C_gene", "IG_LV_gene", "IG_V_gene", "IG_J_gene", "IG_C_gene" , "IG_D_gene")]
#protein_genes <- gene_list$external_gene_name[gene_list$gene_biotype %in% c("protein_coding")]

selected_genes <- rownames(OV4A)[rowSums(OV4A) > 5]
selected_genes <- protein_genes[protein_genes %in% selected_genes]
selected_genes <- selected_genes[!grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", selected_genes)]


OV4A <- SubsetSTData(OV4A, features = selected_genes)




OV4A <- InputFromTable(infotable = infoTable, 
                     platform =  "Visium", minUMICountsPerGene = 5)


img <- Seurat::Read10X_Image(image.dir = 'Desktop/IJC/datasets/IGTP/4A/images/')
Seurat::DefaultAssay(object = img) <- 'Spatial'
rownames(img@coordinates) <- paste0(rownames(img@coordinates), "_1")
img <- img[colnames(x = ST)]
ST[['OC']] <- img



# Load images into a Seurat objectm, mask and plot images
ST <- LoadImages(ST, verbose = TRUE) %>% MaskImages()
ImagePlot(ST)
OC <- RegionNeighbours(OC, id = 11, verbose = T)
OC$nbs_0

ST.FeaturePlot(ST, features = c("IGLC3"), pt.size = 5)


ST <- SubsetSTData(ST, spots = paste0(colnames(OC), "_1"))

OC@tools$Staffli <- ST@tools$Staffli
OC$
FeatureOverlay(ST, features = "clus")
FeatureOverlay(ST, features = "nbs_1")
FeatureOverlay(ST, features = "experiment")

ST$clus <- "0"

x <- mean(ST@reductions$NMF@cell.embeddings[,2]) + 1.5*sd(ST@reductions$NMF@cell.embeddings[,2])
ST$clus[ST@reductions$NMF@cell.embeddings[,2] > x] <- "1"


Idents(ST) <- ST$experiment
ST <- RegionNeighbours(ST, verbose = T, id = "bs_f15")

FeatureOverlay(ST, features = "nbs_bs_f15", pt.size = 3)
