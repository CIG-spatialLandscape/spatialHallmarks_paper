library(BayesSpace)
library(ggplot2)
library(patchwork)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(Seurat)
library(mclust)
library(STutility)
library(cowplot)
library(biomaRt) # need to install this on this workstation
library(sctransform)
library(harmony)
library(RColorBrewer)


OV5A <- readVisium("Desktop/IJC/datasets/IGTP/5A")
colData(OV5A)$sample_name <- "OV5A"
OV5B <- readVisium("/home/malsibai/rnaseq/data/ST/OV_IGTP/sranger_output/Ovario_5B_outs/outs/")
colData(OV5B)$sample_name <- "OV5B"

# combine
OV_combined <- cbind(OV5A, OV5B, deparse.level = 1)

# QC
mito_genes <- rownames(OV_combined)[grep("^MT-", rownames(OV_combined))]

ribo_genes <- rownames(OV_combined)[grep("^RP[Sl]", rownames(OV_combined))]

OV_combined <- addPerCellQC(OV_combined, flatten = T, subsets = list(mt = mito_genes,
                                                                     ribo = ribo_genes))

p <- plot_grid(plotColData(OV_combined, y = "detected", x = "sample_name", colour_by = "sample_name"), plotColData(OV_combined,  y = "total", x = "sample_name", colour_by = "sample_name"),
               plotColData(OV_combined, y = "subsets_mt_percent", x = "sample_name", colour_by = "sample_name"),
               plotColData(OV_combined, y = "subsets_ribo_percent", x = "sample_name", colour_by = "sample_name"), ncol = 2)


jpeg(sprintf("%s/prefilt_QCs.jpg", "./OV5_combined/"), width = 20, height = 10, units = "in", res = 300)
print(p)
dev.off()

# remove mitochondrial and non-coding genes
dim(OV_combined)
OV_combined <- OV_combined[!grepl("^MT-", rownames(OV_combined)), ]
dim(OV_combined)

genes <- rownames(OV_combined)

write.table(genes, file = sprintf("%s/genes.xls", "./"), quote = FALSE, sep = "\t", row.names = FALSE)


mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <-useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
  mart = mart,
  attributes = c("gene_biotype", "external_gene_name"),
  filter ="external_gene_name",
  values = genes$x,
  uniqueRows = TRUE)

annotLookup <- read.delim("~/rnaseq/expression/ST/OV_IGTP/OV_IGTP_OCT/OV_4A/OV_Combined/OV_combined/annotLookup.xls")

protein_coding <- filter(annotLookup, gene_biotype == "protein_coding")

OV_combined <- OV_combined[protein_coding$external_gene_name, ]
dim(OV_combined)
rownames(OV_combined)
colnames(OV_combined)

# Uniquify colnames
colnames(OV_combined) <- uniquifyFeatureNames(OV_combined$spot, OV_combined$sample_name)
OV_combined$spot <- colnames(OV_combined)
colnames(OV_combined@assays@data$counts) <- colnames(OV_combined)

# Normalize with SCt
OV_SCT <- vst(as(counts(OV_combined), "dgCMatrix"), verbosity = 0)$y

OV_combined.SCT <- as(SummarizedExperiment(list(SCT=OV_SCT)), "SingleCellExperiment")
colData(OV_combined.SCT) <- OV_combined@colData

# QC post normalization

OV_combined.SCT <- addPerCellQC(OV_combined.SCT, flatten = T, assay.type = "SCT")

p <- plot_grid(plotColData(OV_combined.SCT, y = "detected", x = "sample_name", colour_by = "sample_name"),
               plotColData(OV_combined.SCT,  y = "total", x = "sample_name", colour_by = "sample_name"),
               ncol = 2)


jpeg(sprintf("%s/postSCT_QCs.jpg", "./OV5_combined/"), width = 10, height = 5, units = "in", res = 300)
print(p)
dev.off()

# preprocess

OV_combined.SCT <- spatialPreprocess(OV_combined.SCT, log.normalize = F, n.PCs = 30, assay.type = "SCT", platform = "Visium")

# Batch correction

OV_combined.SCT <- runUMAP(OV_combined.SCT, dimred = "PCA")
colnames(reducedDim(OV_combined.SCT, "UMAP")) = c("UMAP1", "UMAP2")

p <- ggplot(data.frame(reducedDim(OV_combined.SCT, "UMAP")),
            aes(x = UMAP1, y = UMAP2, color = factor(OV_combined.SCT$sample_name))) +
  geom_point() +
  labs(color = "Sample") +
  theme_bw()

jpeg(sprintf("%s/postSCT_UMAP_preIntegration.jpg", "./OV5_combined/"), width = 10, height = 10, units = "in", res = 300)
print(p)
dev.off()

saveRDS(OV_combined.SCT, "./OV5_combined/OV5_combined.SCT.rds")

### Run Harmony # run it on my PC then get it back here

#OV_combined.SCT <- RunHarmony(OV_combined.SCT, "sample_name", verbose = F)

OV_combined.SCT <- readRDS("./OV5_combined/OV5_combined.SCT.rds")

OV_combined.SCT <- runUMAP(OV_combined.SCT, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(OV_combined.SCT, "UMAP.HARMONY")) = c("UMAP1", "UMAP2")

p <- ggplot(data.frame(reducedDim(OV_combined.SCT, "UMAP.HARMONY")),
            aes(x = UMAP1, y = UMAP2, color = factor(OV_combined.SCT$sample_name))) +
  geom_point() +
  labs(color = "Sample") +
  theme_bw()

jpeg(sprintf("%s/postSCT_UMAP_postIntegration.jpg", "./OV5_combined/"), width = 10, height = 10, units = "in", res = 300)
print(p)
dev.off()

clusterPlot(OV_combined.SCT, "sample_name", color = NA) + #make sure no overlap between samples
  labs(fill = "Sample", title = "Offset check")

OV_combined.SCT$row[OV_combined.SCT$sample_name == "OV5A"] =
  100 + OV_combined.SCT$row[OV_combined.SCT$sample_name == "OV5A"]

clusterPlot(OV_combined.SCT, "sample_name", color = NA) + #make sure no overlap between samples
  labs(fill = "Sample", title = "Offset check")

# find q for clusters
OV_combined.SCT <- qTune(OV_combined.SCT, qs= seq(2,20), platform= "Visium", d = 30)
qPlot(OV_combined.SCT)

q = 8
d = 30

Y <- reducedDim(OV_combined.SCT, "PCA")[, seq_len(d)]
set.seed(101)
init <- Mclust(Y, q, "EEE", verbose = F)$classification

## Run the bayespsace clustering
set.seed(100)
OV_combined.SCT <- spatialCluster(OV_combined.SCT, q=q, d=d, platform = "Visium", init = init,
                                  nrep = 50000, gamma = 3)

OV_combined.SCT$spatial.cluster <- as.factor(OV_combined.SCT$spatial.cluster)

palette <- RColorBrewer::brewer.pal(q, "Set1")

p <- clusterPlot(OV_combined.SCT, palette = palette, color = NA)

pdf(file = "./OV5_combined/BayesSpace_clusters_wo_HE.pdf", width = 10)
print(p)
dev.off()