
#author: Sergi Cervulla
#date: 17/2/2022

#### Packages ####
library(STutility)



### Load data
#create table with the file paths (STutlity format)
infoTable <- data.frame(samples="Desktop/IJC/datasets/IGTP/4A/filtered_feature_bc_matrix.h5",
                        spotfiles="Desktop/IJC/datasets/IGTP/4A/images/tissue_positions_list.csv",
                        imgs="Desktop/IJC/datasets/IGTP/4A/images/tissue_hires_image.png",
                        json="Desktop/IJC/datasets/IGTP/4A/images/scalefactors_json.json")

#create Seurat object with filtering
OC_4A <- InputFromTable(infotable = infoTable, 
                              min.gene.count = 100, 
                              min.gene.spots = 5,
                              min.spot.count = 500,
                              platform =  "Visium")

‘imager’, ‘raster’, ‘magick’, ‘imagerExtra’, ‘ggiraph’, ‘spdep’ are not available for package ‘STutility’



#count matrix
expr.data <- Seurat::Read10X_h5(filename = "Desktop/IJC/datasets/OC_public/filtered_feature_bc_matrix.h5")


OCPublic <- Seurat::CreateSeuratObject(counts = expr.data, project = 'oc_public', assay = 'Spatial', na.rm=T)

img <- Seurat::Read10X_Image(image.dir = 'Desktop/IJC/datasets/OC_public/spatial/')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = OCPublic)]
OCPublic[['OCPublic']] <- img


expr.data <- Seurat::Read10X_h5(filename = "Desktop/IJC/datasets/IGTP/4A/filtered_feature_bc_matrix.h5")
OC4 <- Seurat::CreateSeuratObject(counts = expr.data, project = '4A', assay = 'Spatial', na.rm=T)

img <- Seurat::Read10X_Image(image.dir = 'Desktop/IJC/datasets/IGTP/4A/images/')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = OC4)]
OC4[['OC4A']] <- img



### Filtering
# remove genes with fewer than 10 counts across the whole dataset 
# remove mitochondrial protein coding genes as well as all non-coding RNA biotypes (“antisense,” “lincRNA,” “pseudogenes”)
gene_list <- read.table("Desktop/annotLookup.xls", header = T)
protein_genes <- gene_list$external_gene_name[gene_list$gene_biotype %in% c("protein_coding", "TR_V_gene", "TR_D_gene", "TR_J_gene", "TR_C_gene", "IG_LV_gene", "IG_V_gene", "IG_J_gene", "IG_C_gene" , "IG_D_gene")]


selected_genes <- rownames(OCPublic)[rowSums(OCPublic) > 5]
selected_genes <- protein_genes[protein_genes %in% selected_genes]
selected_genes <- selected_genes[!grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", selected_genes)]
OCPublic <- OCPublic[selected_genes,]  

selected_genes <- rownames(OC4)[rowSums(OC4) > 5]
selected_genes <- protein_genes[protein_genes %in% selected_genes]
selected_genes <- selected_genes[!grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", selected_genes)]
OC4 <- OC4[selected_genes,]  


mito.genes <- grep("^MT-", rownames(expr.data), value = TRUE)


ggplot(a, aes(x=b, y=log(val))) + geom_violin() + geom_jitter(aes(size=0.1)) + scale_size_area(max_size =0.1)


#### SCTtransform
OC.st <- SCTransform(OC.st,return.only.var.genes = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.1, assay = "Spatial")

OC.st <- RunPCA(OC.st, assay = "SCT", verbose = FALSE)

pca=OC.st@reductions$pca
eival= (pca@stdev)^2
varexplaiend = eival/sum(eival)
plot(1:50, cumsum(varexplaiend), type="l")

OC.st <- FindNeighbors(OC.st, reduction = "pca", dims = 1:50)
OC.st <- FindClusters(OC.st, verbose = FALSE, resolution = 0.7)
OC.st <- RunUMAP(OC.st, reduction = "pca", dims = 1:50)

SpatialDimPlot(OC.st, repel=TRUE,label = T)

OC.st_4 <- RunNMF(OC.st, nfactors = 4)
OC.st_4 <- RunNMF(OC, nfactors = 4)

OC.st_4 <- AddMetaData(OC.st_4, as.data.frame(OC.st_4@reductions$NMF@cell.embeddings))

p <- SpatialFeaturePlot(OC.st_4, features = colnames(OC.st_4@meta.data)[6:10], ncol=2)
pdf("Desktop/IJC/datasets/OC_public/figures/factors4_merged2.pdf", width = 5, height = 20)
plot(p)
dev.off()

OC.st_20 <- RunNMF(OC.st, nfactors = 15)
OC.st_20 <- RunNMF(OC, nfactors = 15)

OC.st_20 <- AddMetaData(OC.st_20, as.data.frame(OC.st_20@reductions$NMF@cell.embeddings))

p <- SpatialFeaturePlot(OC.st_20, features = colnames(OC.st_20@meta.data)[6:36], ncol=3)
pdf("Desktop/IJC/datasets/OC_public/figures/factors15/factors15_merged2.pdf", width = 10, height = 60)
plot(p)
dev.off()


FactorGeneLoadingPlot(OC.st_4, factor = "4", topn = 30)  

FactorGeneLoadingPlot(OC.st_20, factor = "14", topn = 30)


################
#Selection top genes
################
library(gprofiler2)

select_top <- function(st, factor) {
  
  vals <- st@reductions$NMF@feature.loadings[, factor]
  logvals <- log(vals)
  logvals <- logvals[is.finite(logvals)]
  thr <- mean(logvals) + 1.645*sd(logvals)
  names(logvals[logvals > thr])
}

i <- 15

top.genes <- list()

for (j in 1:i) {
  top.genes[[j]] <- select_top(OC.st_20, j)
}

top <- select_top(OC.st_20, 5)

a <- gost(query = top, organism = "hsapiens", sources = "GO:BP", evcodes = TRUE)

pathways <- lapply(seq_along(top.genes), function(i) {
  gset <- top.genes[[i]]
  df1 <- gost(query = gset, organism = "hsapiens", sources = "GO:BP", evcodes = TRUE, )$result
  df2 <- gost(query = gset, organism = 'gp__aJd5_EBf2_fFc', evcodes = TRUE)$result
  df <- rbind(df1, df2)
  if (is.null(df)) return(NULL)
  df$factor <- paste0("factor_", i)
  return(df)
})

pathways <- do.call(rbind, pathways)

pathways$factor <- factor(pathways$factor, paste0("factor_", 1:i))
pathways.summarized <- pathways %>% 
  group_by(factor) %>%
  top_n(n = 5, wt = -log10(p_value)) %>%
  mutate(GeneRatio = intersection_size/query_size) %>%
  arrange(factor(factor, levels = paste0("factor_", i:1)), -p_value) %>%
  ungroup() %>%
  mutate(ord = 1:n()) 

sm <- pathways.summarized %>% group_by(factor) %>% summarize(x_start = min(ord)) %>% as.data.frame()
rownames(sm) <- sm$factor
sm <- sm[paste0("factor_", i:1), ]
sm$x_end <- sm$x_start + c(diff(sm$x_start), nrow(pathways.summarized) - max(sm$x_start) + 1)
sm$x_start <- sm$x_start - 0.5
sm$x_end <- sm$x_end - 0.5
sm$x_mid <- (sm$x_end + sm$x_start)/2

# Remove term prefix
pathways.summarized$term_name[pathways.summarized$source == "h.all.v7.1.symbols"] <- pathways.summarized$term_id[pathways.summarized$source == "h.all.v7.1.symbols"]
pathways.summarized$term_name <- gsub(pattern = "_", replacement = " ", x = pathways.summarized$term_name)
pathways.summarized$term_name <- tolower(pathways.summarized$term_name)
pathways.summarized$face <- ifelse(pathways.summarized$source == "GO:BP", "plain", "bold")


p <- ggplot() +
  geom_rect(data = sm, aes(xmin = x_start, xmax = x_end, ymin = 0, ymax = 30, fill = factor), alpha = 0.5, color = "black") +
  geom_segment(data = pathways.summarized, aes(x = ord, xend = ord, y = 0, yend = -log10(p_value)), stat = "identity", linetype = "dashed") +
  geom_point(data = pathways.summarized, aes(ord, -log10(p_value), size = GeneRatio), stat = "identity", shape = 21, fill = "grey") +
  scale_x_continuous(position = "top", breaks = pathways.summarized$ord, labels = pathways.summarized$term_name) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, length.out = 9), expand = c(0, 0)) +
  scale_size_continuous(range = c(0.5, 5)) +
  #guides(fill = FALSE) +
  coord_flip() +
  theme(legend.position = "bottom", 
        axis.text.y = element_text(face = pathways.summarized$face, size = 7), 
        plot.margin = margin(t = 1, r = 1, b = 0, l = 10), 
        axis.text.x = element_text(), 
        panel.background = element_blank(), 
        panel.grid.major.y = element_line(color = "lightgray"), 
        panel.grid.minor.y = element_line(color = "lightgray"), 
        plot.background = element_blank()) +
  labs(fill = "", x = "")
file = paste0("Desktop/IJC/datasets/OC_public/figures/merged15/enrichment.pdf")
pdf(file, width = 12, height = 10)
plot(p)
dev.off()

####
loadings <- OC.st_20@reductions$NMF@feature.loadings[,1]
loadings <- loadings[loadings != 0]
log_loadings <- log(loadings)
log_loadings <- sort(log_loadings, decreasing = T)
smooth <- smth.gaussian(log_loadings)
smooth <- sort(smooth, decreasing = T)
knee  <- uik(1:length(smooth), smooth)


top.genes <- lapply(1:20, function(i) {
  vals <- OC.st_20@reductions$NMF@feature.loadings[, i]
  gg <- data.frame(x = 1:length(vals), vals = sort(vals, decreasing = T))
  gg$gene <- rownames(gg)
  gg <- subset(gg, vals > 0)
  gg$logvals <- log(gg$vals)
  gg$smooth_logvals <- smth.gaussian(x = gg$logvals, windowLength = 10)
  #ggs <- gg[1:500, ]
  ggs <- gg

  knee <- uik(x = ggs$x, y = ggs$smooth_logvals)
  ggs[1:knee, ]$gene
})

l <- c()
for (i in 1:20) {
  l <- c(l,length(top.genes[[i]]))
}

l2 <- c()
for (i in 1:20) {
  l2 <- c(l2,length(top.genes[[i]]))
}

l3 <- c()
for (i in 1:20) {
  l3 <- c(l3,length(top.genes[[i]]))
}

vals <- OC.st_20@reductions$NMF@feature.loadings[, 2]
gg <- data.frame(x = 1:length(vals), vals = sort(vals, decreasing = T))
gg$gene <- rownames(gg)
gg <- subset(gg, vals > 0)
gg$logvals <- log(gg$vals)
gg$smooth_logvals <- smth.gaussian(x = gg$logvals, windowLength = 10)
ggs <- gg[1:500, ]
knee <- uik(x = ggs$x, y = ggs$smooth_logvals)
ggs[1:knee, ]$gene


for (i in 1:20) {
  print(length(top.genes[[i]]))
}


pathways$GeneRatio <- pathways$intersection_size/pathways$query_size
pathways.summarized <- pathways %>% 
  group_by(factor) %>%
  top_n(n = 30, wt = -log10(p_value))
pathways.summarized$term_name[pathways.summarized$source == "h.all.v7.1.symbols"] <- pathways.summarized$term_id[pathways.summarized$source == "h.all.v7.1.symbols"]
p.list <- lapply(1:20, function(i) {
  if (!paste0("factor_", i) %in% pathways.summarized$factor) return(NULL)
  g <- subset(pathways.summarized, factor %in% paste0("factor_", i))
  g$term_name <- ifelse(nchar(g$term_name) > 40, paste0(substr(g$term_name, start = 1, stop = 40), "..."), g$term_name)
  ggplot() +
    geom_point(data = g, aes(reorder(term_name, -log10(p_value)), -log10(p_value), fill = -log10(p_value), size = GeneRatio), color = "black", shape = 21) +
    coord_flip() +
    facet_grid(~factor) +
    scale_size_continuous(range = c(0.5, 8)) +
    scale_fill_gradientn(colours = viridis::magma(n = 9) %>% rev()) +
    theme_minimal() +
    labs(x = "term", y = "")
})
p.list[[10]]

for (i in 1:15){
  file_name = paste0("Desktop/IJC/datasets/OC_public/figures/merged15/pathways/pathways_",i,".jpeg")
  ggsave(file_name, plot = p.list[[i]], height = 7, width = 7)
}

for (i in 1:15) {
  file_name = paste0("Desktop/IJC/datasets/IGTP/5merged/figures/factors15/loadings/loadings_factor", i, ".jpeg")
  ggsave(file_name, plot = FactorGeneLoadingPlot(OC5_20, factor = i, topn = 30))
}



##############

add_path_activities(OC.st_20, assay = "SCT",species = "human", top = 1000)

model <- progeny::getModel(organism = "Human", top = 1000)
common_genes <- intersect(rownames(GetAssayData(OC.st_20, assay = "SCT")), rownames(model))
progeny_scores <- t(model)[, common_genes] %*% GetAssayData(OC.st_20, assay = "SCT")[common_genes, ]
progeny_scores <- scale(t(as.matrix(progeny_scores)))

OC.st_20[['progeny']] <- CreateAssayObject(counts = t(progeny_scores))
DefaultAssay(OC.st_20) <- "progeny"
#pearson correlation
df <- data.frame(row.names = 1:15)

for (factor in 1:15) {
  for (pathway in rownames(OC.st_20)) {
    df[factor,pathway] <-  cor(OC.st_20@reductions$NMF@cell.embeddings[,factor], OC.st_20@assays$progeny@counts[pathway,])
  }
}

myColor = colorRampPalette(c("Darkblue", "white","red"))(100)

pheatmap(as.matrix(df),fontsize=14, 
         fontsize_row = 10, 
         color=myColor,  
         main = "PROGENy (1000)", angle_col = 45,
         treeheight_col = 0,  border_color = NA)

##################################################
samples.combined <- readRDS("Desktop/IJC/datasets/olbrecht2021/RDS/integrated_olalekan&olbrecht.rds")
DefaultAssay(OC.st_20) <- "SCT"
DefaultAssay(samples.combined) <- "SCT"
anchors <- FindTransferAnchors(reference = samples.combined, query = OC.st_20, normalization.method = "SCT", n.trees=1000, reduction = "cca", features = rownames(samples.combined))

predictions.assay <- TransferData(anchorset = anchors, refdata = Idents(samples.combined), prediction.assay = TRUE,
                                  weight.reduction = "cca", dims = 1:30, n.trees = 1000)

OC.st_20[["prediction"]] <- predictions.assay
DefaultAssay(OC.st_20) <- "prediction"
p <- SpatialFeaturePlot(object = OC.st_20, features = rownames(OC.st_20), , ncol = 3)
pdf(file = "Desktop/IJC/datasets/IGTP/4A/figures/seurat_prediction.pdf", width = 15, height = 40)
plot(p)
dev.off()

df <- data.frame(row.names = 1:20)

for (factor in 1:20) {
  for (celltype in rownames(OC.st_20)) {
    df[factor,celltype] <-  cor(OC.st_20@reductions$NMF@cell.embeddings[,factor], OC.st_20@assays$prediction@data[celltype,])
  }
}


df_cell2loc <- data.frame(row.names = 1:20)
for (factor in 1:20) {
  for (celltype in colnames(m)) {
    df_cell2loc[factor,celltype] <-  cor(OC.st_20@reductions$NMF@cell.embeddings[rownames(m),factor], m[,celltype])
  }
}

pheatmap(as.matrix(df_cell2loc),fontsize=14, 
         fontsize_row = 10, 
         color=myColor,  
         main = "Cell2loc", angle_col = 45,
         treeheight_col = 0,  border_color = NA)


df <- data.frame(row.names = rownames(OC.st_20@assays$progeny@counts))

for (pathway in rownames(OC.st_20@assays$progeny@counts)) {
  for (celltype in rownames(OC.st_20@assays$prediction@data)) {
    df[pathway,celltype] <-  cor(OC.st_20@assays$progeny@counts[pathway,], OC.st_20@assays$prediction@data[celltype,])
  }
}


######

p <- SpatialFeaturePlot(OC.st_20, features = SpatiallyVariableFeatures(OC.st_20)[1:50], ncol = 5)
pdf("Desktop/IJC/datasets/IGTP/4A/figures/spatiallyfeatures.pdf", width = 10, height = 20)
plot(p)
dev.off()
##########

library(mistyR)
library(future)

# Seurat
library(Seurat)

# data manipulation
library(Matrix)
library(tibble)
library(dplyr)
library(purrr)

# normalization
library(sctransform)

# resource
library(progeny)

# setup parallel execution
options(future.globals.maxSize = 1024^3)
plan(multisession)


# Define assay for each view
view.assays <- list(
  "main" = "SCT",
  "para.1" = "SCT",
  "para.2" = "SCT",
  "para.3" = "SCT",
  "para.4" = "SCT",
  "para.5" = "SCT",
  "para.6" = "SCT",
  "para.7" = "SCT",
  "para.8" = "SCT",
  "para.9" = "SCT",
  "para.10" = "SCT",
  "para.11" = "SCT",
  "para.12" = "SCT",
  "para.13" = "SCT",
  "para.14" = "SCT",
  "para.15" = "SCT",
  "para.16" = "SCT",
  "para.17" = "SCT",
  "para.18" = "SCT",
  "para.19" = "SCT",
  "para.20" = "SCT"
  )

# Define features for each view
view.features <- list(
  "main" = top.genes[[1]][1:50],
  "para.1" = top.genes[[1]][1:50],
  "para.2" = top.genes[[2]][1:50],
  "para.3" = top.genes[[3]][1:50],
  "para.4" = top.genes[[4]][1:50],
  "para.5" = top.genes[[5]][1:50],
  "para.6" = top.genes[[6]][1:50],
  "para.7" = top.genes[[7]][1:50],
  "para.8" = top.genes[[8]][1:50],
  "para.9" = top.genes[[9]][1:50],
  "para.10" = top.genes[10[]][1:50],
  "para.11" = top.genes[[11]][1:50],
  "para.12" = top.genes[[12]][1:50],
  "para.13" = top.genes[[13]][1:50],
  "para.14" = top.genes[[14]][1:50],
  "para.15" = top.genes[[15]][1:50],
  "para.16" = top.genes[[16]][1:50],
  "para.17" = top.genes[[17]][1:50],
  "para.18" = top.genes[[18]][1:50],
  "para.19" = top.genes[[19]][1:50],
  "para.20" = top.genes[[20]][1:50]
)

# Define spatial context for each view
view.types <- list(
  "main" = "intra",
  "para.1" = "para",
  "para.2" = "para",
  "para.3" = "para",
  "para.4" = "para",
  "para.5" = "para",
  "para.6" = "para",
  "para.7" = "para",
  "para.8" = "para",
  "para.9" = "para",
  "para.10" = "para",
  "para.11" = "para",
  "para.12" = "para",
  "para.13" = "para",
  "para.14" = "para",
  "para.15" = "para",
  "para.16" = "para",
  "para.17" = "para",
  "para.18" = "para",
  "para.19" = "para",
  "para.20" = "para"
)

# Define additional parameters (l in the case of paraview)
view.params <- list(
  "main" = NULL,
  "para.1" = 10,
  "para.2" = 10,
  "para.3" = 10,
  "para.4" = 10,
  "para.5" = 10,
  "para.6" = 10,
  "para.7" = 10,
  "para.8" = 10,
  "para.9" = 10,
  "para.10" = 10,
  "para.11" = 10,
  "para.12" = 10,
  "para.13" = 10,
  "para.14" = 10,
  "para.15" = 10,
  "para.16" = 10,
  "para.17" = 10,
  "para.18" = 10,
  "para.19" = 10,
  "para.20" = 10
)

misty.out <- "Desktop/IJC/datasets/IGTP/4A/misty/"

source("Desktop/IJC/misty_func.r")
misty.results <- run_misty_seurat(
  visium.slide = OC.st_20,
  view.assays = view.assays,
  view.features = view.features,
  view.types = view.types,
  view.params = view.params,
  spot.ids = NULL, # Using the whole slide
  out.alias = misty.out
) %>%
  collect_results()

for (i in 1:20) {
  print(length(top.genes[[i]][1:50]))
}

#######3

df <- data.frame(row.names = 1:15)

for (factor in 1:15) {
  for (celltype in rownames(OC.st)[-22]) {
    df[factor,celltype] <-  cor(OC.st@reductions$NMF@cell.embeddings[,factor], OC.st@assays$prediction@data[celltype,])
  }
}

df <- data.frame(row.names = rownames(OC.st_20@assays$progeny))

for (pathway in rownames(OC.st_20@assays$progeny)) {
  for (celltype in rownames(OC.st_20@assays$prediction)[-22]) {
    df[pathway,celltype] <-  cor(OC.st_20@assays$progeny@data[pathway,], OC.st_20@assays$prediction@data[celltype,])
  }
}
pheatmap(as.matrix(df),fontsize=14, 
         fontsize_row = 10, 
         color=myColor,  
         main = "Progeny (1000)", angle_col = 45,
         treeheight_col = 0,  border_color = NA)


pheatmap(as.matrix(df),fontsize=14, 
         fontsize_row = 10, 
         color=myColor,  
         main = "Progeny (1000)", angle_col = 45,
         treeheight_col = 0,  border_color = NA, breaks = brks)
legend_breaks  = c(-0.8,-0.4,0,0.4, 0.8)

cols <- colorRampPalette(brewer.pal(6,name="PuOr"))(12)
brks <- seq(-0.6,0.6,length.out=100)  
myColor = colorRampPalette(c("Darkblue", "white","red"))(100)
