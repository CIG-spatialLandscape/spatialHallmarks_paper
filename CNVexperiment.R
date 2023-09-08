##################################################
## Project: Cancer Hallmarks
## Script purpose: Establish control gene signatures based on two different perturbations
## Author: Sergi Cervilla* & Mustafa Sibai*
##################################################
library(BayesSpace)
library(tidyverse)
library(infercnv)
library(Seurat)

####################### Step 1: InferCNV state for each sample ##########################
### this script should be executed for each individual sample

sample <- ""
#assign estiamte scores to each sub-spot in the seurat object
gct <- read.table(paste0("", sample, ".gct"), sep = "\t", header = F, skip = 3)
#cluster the estimate scores and sort them from cancer (1) to TME (5)
labels <- kmeans(t(gct[3,3:ncol(gct)]), centers = 5)$cluster

labels_order <- order(sapply(1:5, function(x){
  mean(t(gct[3,3:ncol(gct)])[labels==x])
}))

new_labels <- function(original, order){
  labels <- original
  #Cancer (lowest)
  labels[original==order[1]] <- "1"
  labels[original==order[2]] <- "2"
  labels[original==order[3]] <- "3"
  labels[original==order[4]] <- "4"
  #TME (highest)
  labels[original==order[5]] <- "5"
  return(labels)
}

ordered_labels <- new_labels(labels,labels_order)
#load sub-spot resolution object (as sce object)
sce <- readRDS("")
sce <- sce[,-1]

sce$ordered_label <- ordered_labels

df <- colData(sce)

max_voting <- function(v) {
  v_1 <- cut(as.numeric(v), breaks = c(0,2.5,3.5, 5), labels = c("Cancer", "Buffer","TME"))
  # Count the occurrences of each value
  vote_counts <- table(v_1)
  
  # Find the maximum count(s)
  max_count <- max(vote_counts)
  # Get the values with the maximum count
  max_votes <- names(vote_counts[vote_counts == max_count])
  if (length(max_votes)>1) return("Buffer")
  return(max_votes)
}

max_voting(df$ordered_label)
estimate_clus_spot <- df %>% as.data.frame() %>%
  group_by(spot.idx) %>% summarise(max=max_voting(ordered_label))

#load spot resolution object (as Seurat object)
STobject <- readRDS("")
STobject$clus <- estimate_clus_spot$max

Idents(STobject) <- "clus"
STobject <- subset(STobject, ident=c("TME", "Cancer"))


path_infer <- "~/projects/Hallmarks/inferCNV/"
if (!file.exists(paste0(path_infer, sample, "/"))) dir.create(paste0(path_infer, sample, "/"))

STobject$spot <- colnames(STobject)
write.table(STobject@meta.data[,c("spot","clus"), drop=F], paste0(path_infer, sample, "/annot.tsv"), sep = "\t", row.names = F,
            col.names = F, quote = F)  
#create infercnv object
object_infCNV <- infercnv::CreateInfercnvObject(raw_counts_matrix=as.matrix(STobject@assays$RNA@counts),
                                                gene_order_file=paste0(path_infer, "genes_order.txt"),
                                                annotations_file=paste0(path_infer, sample, "/annot.tsv"),
                                                delim="\t",
                                                ref_group_names=c("TME"),
                                                chr_exclude = c("chrM"))

options(scipen = 100)
rm(STobject, sce, gct)
#run infercnv on spot resolution
object_infCNV = run(object_infCNV,  num_threads = 1,
                    cutoff=0.1, #(see infercnv::run documentation)
                    out_dir=paste0(path_infer, sample, "/output"),
                    cluster_by_groups=F, plot_steps=T, denoise=T,  no_prelim_plot=F, k_obs_groups = 1, HMM=T, leiden_resolution = 0.005) #denoising applies noise reduction for the plot

####################### Step 2: transfer InferCNV from spot to sub-spot resolution ##########################

for (sample in samples)  { 
  print(sample)
  if (!file.exists( paste0("", sample, ""))) {
    cell_groupings <- read.table(paste0('', sample, '/output/infercnv.15_tumor_subclusters.observation_groupings.txt'),header = T)
    cell_groupings$class <- gsub("all_observations_","",cell_groupings$Dendrogram.Group)
    
    cell_groupings$group <- NA
    for (group in names(cnv_grouping[[sample]])) {
      cell_groupings$group[cell_groupings$class %in% cnv_grouping[[sample]][[group]]] <- group
    }
    #load hallmarks scores of sub-spot resoltion
    hallmarks <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hallmarks_scores/", sample, "_hallmarks.txt"))
    hallmarks
    
    
    cnv_clusters <- cell_groupings
    #single cell experiment sub-spot resolution object without imputed genes
    sce <- readRDS(paste0("", sample, "_sce_enhanced.rds"))
    #single cell experiment spot resolution object 
    sce_spot <- readRDS(paste0("Desktop/IJC/datasets/Public/", sample, "/RDS/", sample, ".rds"))
    sce_spot <- SingleCellExperiment(assays =list(SCT = as.matrix(GetAssayData(sce_spot, slot = "data"))), colData = sce_spot@images[[1]]@coordinates)
    meta <- as.data.frame(colData(sce))  
    
    
    for (spot.id in rownames(cnv_clusters)) {
      coord <- colData(sce_spot)[spot.id, c("row", "col")]
      spots <- rownames(meta[meta$spot.row==coord[,1] & meta$spot.col==coord[,2],])
      
      hallmarks[spots, "cnv_cluster"] <- cnv_clusters[spot.id,"group"]
    }
    hallmarks[is.na(hallmarks$cnv_cluster), "cnv_cluster"] <- "Not tumor"
    #save CNV cluster and hallmark score for each sub-spot in a sample
    write.table(hallmarks, paste0("", sample, "_CNV_2.txt"))
    gc()
  }
}

####################### Step 3: Compute genomic distance between pair of clones within a sample ##########################

n_genes <- sapply(samples, function(sample) {
  genes <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/inferCNV/", sample, "/output/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.genes_used.dat"), header = T)
  nrow(genes)
})
for (sample in samples) {
  print(sample)
  #load CNV cluster and hallmark score for each sub-spot in a sample (step 2)
  hallmarks <- read.table(paste0("", sample, "_CNV_2.txt"))
  data <- read.table(paste0("", sample, "/output/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat"), header = T)
  #remove TME clusters
  data <- dplyr::filter(data, grepl('all_observations.all_observations*', cell_group_name))
  groups <- read.table(paste0("", sample, "/output/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings"), header = T)
  #remove TME clusters
  groups <- dplyr::filter(groups, grepl('all_observations.all_observations*', cell_group_name))
  groups$cell_group_name <- str_replace_all(groups$cell_group_name, "all_observations.all_observations_", "")
  data$cell_group_name <- str_replace_all(data$cell_group_name, "all_observations.all_observations_", "")
  #remove clusters with less than 60 sub-spots (10 spots)
  tbl <- table(groups$cell_group_name)
  low_cnv <- names(tbl)[table(groups$cell_group_name) < 10]
  groups[groups$cell_group_name %in% low_cnv, "cell_group_name"] <- "excluded"
  hallmarks[hallmarks$cnv_cluster %in% low_cnv, "cnv_cluster"] <- "excluded"
  data[data$cell_group_name %in% low_cnv, "cell_group_name"] <- "excluded"
  
  df <- data.frame(matrix(0,nrow = length(unique(groups$cell_group_name)), ncol = length(unique(data$gene))))
  colnames(df) <- unique(data$gene)
  rownames(df) <- unique(groups$cell_group_name)
  for (row in 1:nrow(data)) {
    grouped_state <- ifelse(data[row, "state"] > 3, 1, -1)
    df[data[row, "cell_group_name"], data[row, "gene"]] <- grouped_state
  }

  df_var <- as.data.frame(matrix(ncol = 16, nrow = 0))
  labels <- rownames(df)
  #compute jaccard index across all pair of clones
  for (i in 1:(length(labels)-1)) {
    for (j in (i+1):length(labels)) {
      
      intersect <- sapply(1:ncol(df), function(x) {
        df[labels[i],x] == df[labels[j],x]
      }) %>% sum() +  n_genes[sample]-length(unique(data$gene))
      h_diff <- sapply(paste0("H", 1:13), function(h) {
        mean(hallmarks[hallmarks$cnv_cluster == labels[i], h], na.rm=TRUE) - mean(hallmarks[hallmarks$cnv_cluster == labels[j], h], na.rm=TRUE)
      })
      df_var <- rbind(df_var, c(h_diff, intersect/n_genes[sample], paste0(labels[i], "_", labels[j]), sample))
    }
  }
  df_var[,1:14] <- apply(df_var[,1:14], 2, as.numeric)
  colnames(df_var) <- c(paste0("H", 1:13), "var", "group", "sample")
  
  write.table(df_var, paste0("", sample, "_jaccard.txt"))
}
