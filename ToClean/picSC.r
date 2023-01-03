samples.combined <- readRDS("Desktop/IJC/datasets/Public/SC/olbrecht2021/RDS/integrated_olalekan&olbrecht.rds")


OV <- Seurat::CreateSeuratObject(Seurat::GetAssayData(samples.combined, slot = "counts"), meta.data = samples.combined@meta.data)

OV@reductions <- samples.combined@reductions
OV@assays$RNA <- samples.combined@assays$RNA 
OV@assays$integrated <- samples.combined@assays$integrated
OV@assays$SCT <- samples.combined@assays$SCT



DimPlot(OV, group.by = "orig.ident", cols = DiscretePalette(13,palette = "alphabet")) + labs(title = "Patients")
DimPlot(OV, group.by = "sample_site", cols = RColorBrewer::brewer.pal(3, name = "Paired")) + labs(title = "Sample site")

OV$sample_site
DiscretePalette(13)
OV$orig.ident[!is.na(OV$patient_id)] <- paste0(OV$patient_id[!is.na(OV$patient_id)], "olbrecht")
OV$sample_site[is.na(OV$patient_id)] <- "omentum"
Idents(OV) <- Idents(samples.combined)

DimPlot(OV, label = T, cols = DiscretePalette(21,palette = "alphabet2")) + labs(title = "Cell type")



DotPlot(OV, assay = "integrated", features = c("CD1C", "CD1E", "CD83", "CD74", "AIF1", "LST1", "C1QB", "CD163", "APOE", "CCL5"))

DotPlot(OV, features = c("C1QB", "CD14", "SPP1"))


DotPlot(OV, assay = "integrated", features = c("CD74", "HLA-DQA1", "CD83",
                                               "CD1E", "AIF1", "LST1",
                                               "C1QB", "CD163", "APOE", 
                                               "CCL5", "CCL4","GZMA",
                                               "PRF1", "GNLY", "KLRD1",
                                               "MKI67", "STMN1", "TOP2A",
                                               "IL7R", "CCR7", "FYB1",
                                               "FOXP3", "IL2RA",  "CTLA4",
                                               "GATA2", "KIT", "MS4A2",
                                               "IRF8", "IRF4", "BCL11A",
                                               "MS4A1", "CD79A", "CD79B",
                                               "MCAM", "ACTA2", "MYH11",
                                               "PECAM1", "VWF", "PALMD",
                                               "IGKC", "IGHG1", "IGHM",
                                               "C7", "MGP", "CXCL14",
                                               "C3", "SERPINE1", "S100A10",
                                               "PDPN", "MMP11", "VCAN",
                                               "WFDC2", "MUC16", "EPCAM"
                                               ), cols = "RdBu", col.min = -3) + theme(axis.text.x = element_text(angle = 90))+ 
  scale_color_viridis_c()
