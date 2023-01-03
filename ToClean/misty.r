# MISTy
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

setwd("Desktop/IJC/datasets/")

# Load the HDF5 object and normalize the expression
seurat.vs <-
  Load10X_Spatial(
    data.dir = "FFPE_PC_IF/",
    filename = "Visium_FFPE_Human_Prostate_IF_filtered_feature_bc_matrix.h5"
  )

sct.data <- vst(GetAssayData(
  object = seurat.vs,
  slot = "counts",
  assay = "Spatial"
),
verbosity = 0
)

seurat.vs[["SCT"]] <- CreateAssayObject(data = sct.data$y)

gene.expression <- GetAssayData(seurat.vs, assay = "SCT")
coverage <- rowSums(gene.expression > 0) / ncol(gene.expression)
slide.markers <- names(which(coverage >= 0.01))

foot.prints <- getModel(top = 15)

#get gens involved in each pathway of progeny
androgen.footprints <- getModel(top = 15) %>%
  rownames_to_column("gene") %>%
  filter(PI3K != 0, gene %in% slide.markers) %>%
  pull(gene)
androgen1.footprints <- getModel(top = 15) %>%
  rownames_to_column("gene") %>%
  filter(Androgen != 0, gene %in% slide.markers) %>%
  pull(gene)

hypoxia.footprints <- getModel(top = 15) %>%
  rownames_to_column("gene") %>%
  filter(Hypoxia != 0, gene %in% slide.markers) %>%
  pull(gene)
tgfb.footprints <- getModel(top = 15) %>%
  rownames_to_column("gene") %>%
  filter(TGFb != 0, gene %in% slide.markers) %>%
  pull(gene)


#define the genes involved in each view
first.footprints <- m[1:30]
second.footprints <- pc_ligands


# Define assay for each view
view.assays <- list(
  "main" = "SCT",
  "para.markers" = "SCT",
  "para.lig" = "SCT",
  "para.rec" = "SCT",
  "para.svg" = "SCT"
)

# Define features for each view
view.features <- list(
  "main" = rownames(markers),
  "para.markers" = rownames(markers),
  "para.lig" = pc_ligands,
  "para.rec" = pc_receptor,
  "para.svg" = sp_features
)

# Define spatial context for each view
view.types <- list(
  "main" = "intra",
  "para.markers" = "para",
  "para.lig" = "para",
  "para.rec" = "para",
  "para.svg" = "para"
)

# Define additional parameters (l in the case of paraview)
view.params <- list(
  "main" = NULL,
  "para.markers" = 10,
  "para.lig" = 10,
  "para.rec" = 10,
  "para.svg" = 10
)
#folder to store results
misty.out <- "misty"

#run misty model
misty.results <- run_misty_seurat(
  visium.slide = seurat.vs,
  view.assays = view.assays,
  view.features = view.features,
  view.types = view.types,
  view.params = view.params,
  spot.ids = NULL, # Using the whole slide
  out.alias = misty.out
) %>%
  collect_results()



#plot results
misty.results %>%
  plot_improvement_stats("gain.R2") %>%
  plot_improvement_stats("gain.RMSE")

misty.results %>% plot_view_contributions()

misty.results %>% plot_interaction_heatmap(view = "intra", 0.5)
misty.results %>% plot_interaction_heatmap(view = "para.1_10", clean = T)
misty.results %>% plot_interaction_heatmap(view = "para.2_10",clean = T, 4)

misty.results %>% plot_contrast_heatmap(viewmisty.results %>% plot_interaction_communities("para.2_10")

                                        plot_interaction_communities(misty.results = MISTy_out, view = "b")

predictors <- list()
for (i in names(misty.results$importances.aggregated)) {
  predictors[[i]] = tidyr::gather(misty.results$importances.aggregated[[i]], 
                             "Predicted",
                             "Importance", -Predictor)
  predictors[[i]] <- predictors[[i]][order(predictors[[i]]$Importance, decreasing = T),]
}

predictors[["intra"]][predictors[["intra"]]$Predicted=="CCN1",] #6
predictors[["para.markers_10"]][predictors[["para.markers_10"]]$Predicted=="CCN1",] #5
predictors[["para.rec_10"]][predictors[["para.rec_10"]]$Predicted=="CCN1",] #PTRPF
predictors[["para.lig_10"]][predictors[["para.lig_10"]]$Predicted=="CCN1",] #8
predictors[["para.svg_10"]][predictors[["para.svg_10"]]$Predicted=="CCN1",] #

features <- predictors[["intra"]][predictors[["intra"]]$Predicted=="CCN1",]$Predictor[1:6]
features[7] <- "CCN1"
SpatialFeaturePlot(prostate_if, features)

features <- predictors[["para.markers_10"]][predictors[["para.markers_10"]]$Predicted=="CCN1",]$Predictor[1:5]
features[6] <- "CCN1"
SpatialFeaturePlot(prostate_if, features)

SpatialFeaturePlot(prostate_if, c("PTPRF", "CCN1"))

features <- predictors[["para.lig_10"]][predictors[["para.lig_10"]]$Predicted=="CCN1",]$Predictor[1:8]
features[9] <- "CCN1"
SpatialFeaturePlot(prostate_if, features, ncol=5)

features <- predictors[["para.svg_10"]][predictors[["para.svg_10"]]$Predicted=="CCN1",]$Predictor[1:9]
features[10] <- "CCN1"
SpatialFeaturePlot(prostate_if, features, ncol=5)




predictors[["intra"]][predictors[["intra"]]$Predicted=="A2M",] #8
predictors[["para.markers_10"]][predictors[["para.markers_10"]]$Predicted=="A2M",] #5
predictors[["para.rec_10"]][predictors[["para.rec_10"]]$Predicted=="A2M",] #PTRPF
predictors[["para.lig_10"]][predictors[["para.lig_10"]]$Predicted=="A2M",] #8
predictors[["para.svg_10"]][predictors[["para.svg_10"]]$Predicted=="A2M",] #

features <- predictors[["intra"]][predictors[["intra"]]$Predicted=="A2M",]$Predictor[1:8]
features[9] <- "A2M"
SpatialFeaturePlot(prostate_if, features)

features <- predictors[["para.markers_10"]][predictors[["para.markers_10"]]$Predicted=="A2M",]$Predictor[1:8]
features[9] <- "A2M"
SpatialFeaturePlot(prostate_if, features)

SpatialFeaturePlot(prostate_if, c("PTPRF", "A2M"))

features <- predictors[["para.lig_10"]][predictors[["para.lig_10"]]$Predicted=="A2M",]$Predictor[1:6]
features[7] <- "A2M"
SpatialFeaturePlot(prostate_if, features, ncol=5)

features <- predictors[["para.svg_10"]][predictors[["para.svg_10"]]$Predicted=="A2M",]$Predictor[1:9]
features[10] <- "A2M"
SpatialFeaturePlot(prostate_if, features, ncol=5)
################################################################################




prostate_if <- FindSpatiallyVariableFeatures(prostate_if, assay = "SCT", features = VariableFeatures(prostate_if), selection.method = "markvariogram")
sp_features <- SpatiallyVariableFeatures(prostate_if, selection.method = "markvariogram")


#################

progeny_if %>%
  plot_improvement_stats("gain.R2") %>%
  plot_improvement_stats("gain.RMSE")

progeny_pc %>%
  plot_improvement_stats("gain.R2") %>%
  plot_improvement_stats("gain.RMSE")

progeny_pc %>% plot_view_contributions()
progeny_if %>% plot_view_contributions()



plot_contrast_results(progeny_if, progeny_pc, views = "intra")
plot_contrast_results(progeny_pc, progeny_if, views = "intra")

plot_interaction_communities(misty.results = progeny_if, view = "intra")
plot_interaction_communities(misty.results = progeny_pc, view = "intra")



MISTy_out_if %>%
  plot_improvement_stats("gain.R2") %>%
  plot_improvement_stats("gain.RMSE")

MISTy_out_pc %>%
  plot_improvement_stats("gain.R2") %>%
  plot_improvement_stats("gain.RMSE")

MISTy_out_if %>% plot_interaction_heatmap(view = "intra", clean = T)
MISTy_out_if %>% plot_interaction_heatmap(view = "para",clean = T, 4)

MISTy_out_pc %>% plot_interaction_heatmap(view = "intra", clean = T)
MISTy_out_pc %>% plot_interaction_heatmap(view = "para",clean = T, 4)

plot_contrast_results(progeny_if, MISTy_out_pc, views = "intra")
plot_contrast_results(progeny_pc, MISTy_out_if, views = "intra")

plot_interaction_communities(misty.results = MISTy_out_pc, view = "intra")
plot_interaction_communities(misty.results = progeny_pc, view = "intra")

MISTy_out_if$importances.aggregated %>% arrange(Predictor)

predictors_if = tidyr::gather(MISTy_out_if$importances.aggregated$para, 
                              "Predicted",
                              "Importance", -Predictor)
predictors_if <- predictors_if[order(predictors_if$Importance, decreasing = T),]
predictors_pc = tidyr::gather(MISTy_out_pc$importances.aggregated$para, 
                              "Predicted",
                              "Importance", -Predictor)
predictors_pc <- predictors_pc[order(predictors_pc$Importance, decreasing = T),]

predictors_if[predictors_if$Predicted == "TGM4" & predictors_if$Importance > 1, ]
predictors_pc[predictors_pc$Predicted == "TGM4" & predictors_if$Importance > 1, ]






#######################################################################