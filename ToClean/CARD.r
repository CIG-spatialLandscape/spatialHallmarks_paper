################################################################################

#                                 Deconvolution

################################################################################

library(CARD)

sc <- readRDS("Desktop/IJC/datasets/Public/SC/olbrecht2021/RDS/integrated_olalekan&olbrecht.rds")
sc_count <- GetAssayData(sc, assay = "RNA", slot = "counts")
batch_var <- "patient_id"
celltype_var <- "MajorTypes"
select_celltypes <- c("CAFs", "TAMs", "EC", "DC", "NK", "Myofibroblast", "Epithelial", "Epithelial_Proliferative", "pDC",
                      "Bcells", "PlasmaB", "Mast_cells", "CD4+Tcells", "CD8+Tcells", "Treg")
paste0("'",unique(sc$MajorTypes), "'", collapse = ", ")
select_celltypes <- c('CAFs', 'EC', 'Epithelial', 'TAMs', 'M1 Macrophages',
'Epithelial_Proliferative', 'Myofibroblast', 'DC', 'PlasmaB', 'NK', 'pDC', 'Bcells', 'Mast_cells')
select_celltypes <- as.character(unique(sc_meta[,2])[c(2,3,5,6,8,9,12,13,14,15,16,18,19,20)])
sc_meta <- sc@meta.data[, c(batch_var, celltype_var)]
rm(sc)
gc()

STobject <- readRDS("Desktop/IJC/datasets/IGTP/4A/RDS/OV4A.rds")
spatial_count <- GetAssayData(STobject, assay = "RNA", slot = "counts")
spatial_location <- STobject@images[[1]]@coordinates[,c("col", "row")]
colnames(spatial_location) <- c("x", "y")
spatial_location[,2] <- -spatial_location[,2]



CARD_obj = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = celltype_var,
  ct.select = select_celltypes,
  sample.varname = batch_var,
  minCountGene = 0,
  minCountSpot = 0) 


CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
CARD_obj@Proportion_CARD
STobject[["Deconvolution"]] <- CreateAssayObject(t(CARD_obj@Proportion_CARD))

CARD.visualize.Cor(CARD_obj@Proportion_CARD,colors = NULL)




coords <- STobject@images[[1]]@coordinates
a <- c()
for (i in colnames(STobject)) {
  x <- sapply(get_neighbors(i), function(spot) {
    as.numeric(STobject@assays$Deconvolution["DC", spot])
  })
  a <- c(a, mean(x))
}

plot(as.numeric(STobject@assays$Deconvolution["Epithelial", ]), a)


as.numeric(STobject@assays$Deconvolution["Tcells", ])






