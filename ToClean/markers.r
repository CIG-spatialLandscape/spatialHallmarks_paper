###########################################################################################

#             Canonical markers tumour cells

# ref: https://www.biorxiv.org/content/10.1101/2020.10.26.354829v1.supplementary-material

###########################################################################################

DefaultAssay(samples.combined) <- "SCT"

## B cells
bcells <- c("MS4A1", "CD79A", "CD83", "CD79B", "CD37", "CD19")
FeaturePlot(samples.combined, features = bcells)
DotPlot(samples.combined, features = bcells)

## Proliferative B cells
prolbcells <- c("CD79B", "MS4A1", "CD79A", "STMN1", "KIAA0101", "TUBA1B", "CD19")
FeaturePlot(samples.combined, features = prolbcells)
DotPlot(samples.combined, features = prolbcells)

## Plasma B cells
plasmab <- c("IGKC", "IGHG1", "IGHM", "TNFRSF17", "SDC1", "CD38")
FeaturePlot(samples.combined, features = plasmab)
DotPlot(samples.combined, features = plasmab)

## Naive T cells
naivetcells <- c("IL7R", "TCF7", "CCR7", "LEF1", "SELL")
FeaturePlot(samples.combined, features = naivetcells)
DotPlot(samples.combined, features = naivetcells)

## Regulatory T cells
regulatoryT <- c("FOXP3", "TNFRSF18", "IL2RA", "TIGIT", "CTLA4", "IKZF2")
FeaturePlot(samples.combined, features = regulatoryT)
DotPlot(samples.combined, features = regulatoryT)

## T helper
thelper <- c("CXCL13", "TNFRSF4", "TNFRSF18", "BATF", "TIGIT", "SOX4", "TNFRSF25", "CTLA4", "RORA", "XCL1", "TNFSF8", "PPIA", "STAT5A", "TOX", "PDCD1")
FeaturePlot(samples.combined, features = thelper)
DotPlot(samples.combined, features = thelper)

## Th17 cells
th17 <- c("TNFSF8", "IL7R", "CXCR4", "VPS37B", "GZMK", "STAT4", "TGFB1", "XCL1", "XCL2", "CCR7", "CRTAM")
FeaturePlot(samples.combined, features = th17)
DotPlot(samples.combined, features = th17)

## Proliferative T cells
proltcells <- c("STMN1", "MKI67", "CDK1")
FeaturePlot(samples.combined, features = proltcells)
DotPlot(samples.combined, features = proltcells)

## Recently activated CD4 T cells
racd4tcells <- c("CCL4", "IFITM1", "CD69", "PRF1", "BCL3", "IL7R", "IFITM3", "TCF7", "CD81", "CXCR4", "GZMK", "GZMM", "IFITM2") #Recently activated CD4 T cells
FeaturePlot(samples.combined, features = racd4tcells)
DotPlot(samples.combined, features = racd4tcells)

## Naive memory CD4 T cells
nmcd4tcells <- c("IL7R", "TCF7", "CCL5", "IFITM1") #Naive memory CD4 T cells
FeaturePlot(samples.combined, features = nmcd4tcells)
DotPlot(samples.combined, features = nmcd4tcells)

## Transitional Memory CD4 T cells
tmcd4tcells <- c("CXCL13", "TNFRSF4", "TIGIT", "IL6ST", "PASK", "KLRB1", "CD40LG","TOX", "LEF1", "ICOS", "CD28", "TOX2", "CCR7", "CD247", "RORA", "PDCD1") #Transitional Memory CD4 T cells
FeaturePlot(samples.combined, features = tmcd4tcells)
DotPlot(samples.combined, features = tmcd4tcells)

## Pre-exhausted CD8 T cells
pecd8tcells <- c("ISG15", "IFI44L", "IFI6", "IFIT3", "IFIT1", "IFI44", "IFI35", "IRF7", "IFIT2", "LAG3", "IFITM1", "IFI16", "IFI27", "IFNG", "GZMB", "GZMK", "PRF1", "HAVCR2", "IFIH1", "GZMA", "IRF9", "CXCL13", "GZMH", "IFIT5", "PDCD1") #Pre-exhausted CD8 T cells
FeaturePlot(samples.combined, features = pecd8tcells)
DotPlot(samples.combined, features = pecd8tcells)

## Cytotoxic CD8 T cells
cytcd8tcells <- c("CCL5", "GZMA", "GZMK", "NKG7", "GNLY", "GZMH", "GZMM", "PRF1", "CXCR3", "GZMB", "CCL4", "IFNG") #Cytotoxic CD8 T cells
FeaturePlot(samples.combined, features = cytcd8tcells)
DotPlot(samples.combined, features = cytcd8tcells)

## Effector memory CD8 T cells
emcd8tcells <- c("GZMM", "IFITM1", "GZMK", "IFNG", "CCL5") #Effector memory CD8 T cells 
FeaturePlot(samples.combined, features = emcd8tcells)
DotPlot(samples.combined, features = emcd8tcells)

## Terminally exhausted CD8 T cells
tecd8tcells <- c("CXCL13", "LAG3", "GZMB", "CCL5", "NKG7", "IFNG", "GZMA", "HAVCR2", "GNLY", "PDCD1", "TIGIT", "TNFRSF9", "ENTPD1", "CTLA4", "PRF1", "TOX", "GZMH", "GZMK") #Terminally exhausted CD8 T cells
FeaturePlot(samples.combined, features = tecd8tcells)
DotPlot(samples.combined, features = tecd8tcells)

## NK
nk <- c("PRF1", "GNLY", "KLRD1", "KLRF1", "GZMH", "GZMB", "KLRB1", "GZMA", "GZMM", "CD160", "CD244", "KLRC1", "NCR1")
FeaturePlot(samples.combined, features = nk)
DotPlot(samples.combined, features = nk)

## SPP1 TAMs
spp1tam <- c("SPP1", "APOE", "SEPP1", "MMP9", "CD163")
FeaturePlot(samples.combined, features = spp1tam)
DotPlot(samples.combined, features = spp1tam)

## M2 TAMs
m2tam <- c("C1QB", "APOE", "C1QA", "C1QC", "APOC1", "SEPP1", "SPP1", "CD163")
FeaturePlot(samples.combined, features = m2tam)
DotPlot(samples.combined, features = m2tam)

## proinflomatory TAMs
protam <- c("CXCL8", "IL1B", "S100A9", "S100A8", "CCL2", "IL8", "CD68", "IL6", "IL1A") #proinflomatory TAMs
FeaturePlot(samples.combined, features = protam)
DotPlot(samples.combined, features = protam)

## Proliferative monocytes and macrophages
pmonmac <- c("C1QB", "C1QC", "SPP1", "C1QA", "TUBA1B", "STMN1", "APOE", "APOC1", "MKI67", "CD14", "KIAA0101", "TOP2A", "CD68", "SPI1", "BIRC5", "CSF1R", "PCNA", "FCER1A", "IL1B", "CD1C")
FeaturePlot(samples.combined, features = pmonmac)
DotPlot(samples.combined, features = pmonmac)

## Monocytes
monocytes <- c("S100A8", "S100A9", "FCN1", "VCAN", "AIF1", "SPI1", "CD14", "APOBEC3A", "CSF1R", "ASAH1")
FeaturePlot(samples.combined, features = monocytes)
DotPlot(samples.combined, features = monocytes)

## cDC
cDC <- c("FSCN1", "HLA.DRA", "HLA.DRB1", "SPI1", "CLEC9A", "XCR1", "ITGAX")
FeaturePlot(samples.combined, features = cDC)
DotPlot(samples.combined, features = cDC)

## pDC
pDC <- c("TCF4", "IRF8", "IRF4", "BCL11A", "SPIB", "CLEC4C", "RUNX2")
FeaturePlot(samples.combined, features = pDC)
DotPlot(samples.combined, features = pDC)

## mDC
mDC <- c("HLA.DRA", "HLA.DRB1", "SPI1", "CD68", "CD83", "ITGAX", "CD1D")
FeaturePlot(samples.combined, features = mDC)
DotPlot(samples.combined, features = mDC)

## Mast cells
mastcells <- c("TPSB2", "TPSAB1", "CPA3", "GATA2", "KIT", "MS4A2")
FeaturePlot(samples.combined, features = mastcells)
DotPlot(samples.combined, features = mastcells)


markers <- list()
markers[["B cells"]] <- bcells
markers[["Proliferative B cells"]] <- prolbcells
markers[["Plasma B cells"]] <- plasmab
markers[["Naive T cells"]] <- naivetcells
markers[["Regulatory T cells"]] <- regulatoryT 
markers[["T helper"]] <- thelper
markers[["Th17 cells"]] <- th17
markers[["Proliferative T cells"]] <- proltcells
markers[["Recently activated CD4 T cells"]] <- racd4tcells
markers[["Naive memory CD4 T cells"]] <- naivetcells
markers[["CD4 transitional memory "]] <- tmcd4tcells
markers[["CD8 pre-exhausted"]] <- pecd8tcells
markers[["CD8  Cytotoxic T cells"]] <- cytcd8tcells
markers[["CD8 effector memory "]] <-  emcd8tcells
markers[["CD8 terminally exhausted"]] <- tecd8tcells
markers[["NK"]] <- nk
markers[["Macrophages spp1"]] <- spp1tam
markers[["TAMs C1QC"]] <- m2tam
markers[["TAMs proinflamatory"]] <- protam
markers[["Macro. and mono. prolif"]] <- pmonmac
markers[["Monocytes"]] <- monocytes
markers[["cDC"]] <- cDC
markers[["pDC"]] <- pDC
markers[["mDC"]] <- mDC
markers[["Mast cells"]] <- mastcells
 
