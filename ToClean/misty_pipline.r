library(OmnipathR)
library(tidyverse)
library(Seurat)
library(cowplot)
library(mistyR)
library(future)

#load functions for misty
setwd("Desktop/IJC/datasets/")
source("../misty_func.r")

plan(multisession)

#read seurat object containing spatial data
prostate_if <- readRDS("pc_seurat.rds")

#get ligands and receptors from omniPath
lig_rec = import_intercell_network(
  interactions_param = list(datasets = c('ligrecextra', 'omnipath', 'pathwayextra'),
                            transmitter_param = list(parent = 'ligand'),
                            receiver_param = list(parent = 'receptor')
  ))

ligands = unique(lig_rec$source_genesymbol)


pc_ligands = ligands[ligands %in% rownames(prostate_if)] 
gex = as.matrix(prostate_if@assays$SCT@data)[pc_ligands,]
lig_coverage = rowSums(gex>0)/ncol(gex)
pc_ligands = names(lig_coverage[lig_coverage>=.3])

#create an assay storing the score of each pathway in ST data
progeny_scores = progeny::progeny(expr = as.matrix(prostate_if[["SCT"]]@data),
                                  scale=TRUE, 
                                  organism="Human", top=1000, perm=1)

prostate_if[['progeny']] = CreateAssayObject(counts = t(progeny_scores))

#values for the paramater l (radius which defines the paraview) / optimal will be selected
ls = c(2,5,10,20,50,100)

misty_out_path = "misty" # Defining the output folder alias

#run misty model
test_opath = lapply(ls,para_ppln_seurat,
                    visium_slide = prostate_if,
                    intra_assay = "SCT",
                    intra_features = NULL,
                    para_assay = "SCT",
                    para_features = sp_features,
                    spot_ids = NULL,
                    out_alias = misty_out_path)

#loading results / aggregation 
get_optimal(out_dir_name = "misty",ls = ls)
MISTy_out = MISTy_aggregator(results_folder = "misty_optim")

#plotting results
plot_misty_performance(MISTy_out = MISTy_out)

plot_misty_importance(MISTy_out, portance_cut = 1)

#check list of predictors
tbl <- plot_misty_importance(MISTy_out, importance_cut = 1)
tbl[[2]][tbl[[2]]$Predicted == "TGM4",]$Predictor[1:10] #check list of predictors
SpatialFeaturePlot(prostate_if, feature =tbl[[2]][tbl[[2]]$Predicted == "TGFb",]$Predictor[1:6] )

#test
MISTy_out$importances.aggregated <- MISTy_out$importance
MISTy_out$importances.aggregated[["b"]] <- MISTy_out$importances.aggregated[[1]]
MISTy_out$importances.aggregated[["a"]] <- MISTy_out$importances.aggregated[[2]]
plot_interaction_communities(misty.results = MISTy_out, view = "b")