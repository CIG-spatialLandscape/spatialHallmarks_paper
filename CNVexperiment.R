library(stringr)
library(tidyverse)


##################################################
## Project: Cancer Hallmarks
## Script purpose: Compute Radar scores of TME Hallmarks within Cancer spots
## Author: Sergi Cervilla* & Mustafa Sibai*
##################################################

samples <- list.files("")

## This loop will assign CNV clusters from spot-resolution to sub-spot resolution
for (sample in samples)  { 
  print(sample)
  if (T) {
    cell_groupings <- read.table(paste0('', sample, '/output_bayes/infercnv.19_HMM_pred.Bayes_Net.Pnorm_0.2.observation_groupings.txt'),header = T)
    cell_groupings$class <- gsub("all_observations_","",cell_groupings$Dendrogram.Group)
    
    #load hallmarks scores of sub-spot resoltion
    hallmarks <- read.table(paste0("/df/", sample, "_hallmarks.txt"))
    
    cnv_clusters <- cell_groupings
    #single cell experiment sub-spot resolution object without imputed genes
    sce <- readRDS(paste0("/sce/", sample, "_sce_enhanced.rds"))
    #single cell experiment spot resolution object 
    sce_spot <- readRDS(paste0("/RDS/", sample, ".rds"))
    sce_spot <- SingleCellExperiment(assays =list(SCT = as.matrix(GetAssayData(sce_spot, slot = "data"))), colData = sce_spot@images[[1]]@coordinates)
    meta <- as.data.frame(colData(sce))  
    
    for (spot.id in rownames(cnv_clusters)) {
      coord <- colData(sce_spot)[spot.id, c("row", "col")]
      spots <- rownames(meta[meta$spot.row==coord[,1] & meta$spot.col==coord[,2],])
      
      hallmarks[spots, "cnv_cluster"] <- cnv_clusters[spot.id,"class"]
    }
    hallmarks[is.na(hallmarks$cnv_cluster), "cnv_cluster"] <- "Not tumor"
    #save CNV cluster and hallmark score for each sub-spot in a sample
    write.table(hallmarks, paste0("/df/", sample, "_CNV_bayes.txt"))
  }
}

# Number of genes used in each sample
n_genes <- sapply(samples, function(sample) {
  genes <- read.table(paste0("", sample, "/output_bayes/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.2.pred_cnv_genes.dat"), header = T)
  nrow(genes)
})



## This loop will compute the genomic distance between clones from the same sample
## This loop will collapse pair of clones that are at least 99% similar in order to avoid the overclustering coming from the leiden resolution
for (sample in samples) {
  print(sample)
  #load CNV cluster and hallmark score for each sub-spot in a sample (step 2)
  hallmarks <- read.table(paste0("/df/", sample, "_CNV_bayes.txt"))
  data <- read.table(paste0("", sample, "/output_bayes/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.2.pred_cnv_genes.dat"), header = T)
  #remove TME clusters
  data <- dplyr::filter(data, grepl('all_observations.all_observations*', cell_group_name))
  groups <- read.table(paste0("", sample, "/output_bayes/infercnv.19_HMM_pred.Bayes_Net.Pnorm_0.2.observation_groupings.txt"), header = T)
  #remove TME clusters
  groups <- dplyr::filter(groups, grepl('all_observations_*', Dendrogram.Group))
  groups$cell_group_name <- str_replace_all(groups$Dendrogram.Group, "all_observations_", "")
  data$cell_group_name <- str_replace_all(data$cell_group_name, "all_observations.all_observations_", "")
  #remove clusters with less than 60 sub-spots (10 spots)
  tbl <- table(groups$cell_group_name)
  low_cnv <- names(tbl)[table(groups$cell_group_name) < 10]
  groups[groups$cell_group_name %in% low_cnv, "cell_group_name"] <- "excluded"
  hallmarks[hallmarks$cnv_cluster %in% low_cnv, "cnv_cluster"] <- "excluded"
  data[data$cell_group_name %in% low_cnv, "cell_group_name"] <- "excluded"
  groups <- groups %>% filter(cell_group_name != "excluded")
  hallmarks <- hallmarks %>% filter(cnv_cluster != "excluded")
  data <- data %>% filter(cell_group_name != "excluded")
  
  #get the CNV states for each clone
  df <- data.frame(matrix(0,nrow = length(unique(groups$cell_group_name)), ncol = length(unique(data$gene))))
  colnames(df) <- unique(data$gene)
  rownames(df) <- unique(groups$cell_group_name)
  for (row in 1:nrow(data)) {
    grouped_state <- ifelse(data[row, "state"] > 3, 1, -1)
    df[data[row, "cell_group_name"], data[row, "gene"]] <- grouped_state
  }
  
  #identify the clusters that have exactly the same CNV profile
  df_uq <- unique(df)
  duplicated.list <- df %>% group_by(across(everything())) %>% rownames_to_column(var = "clone") %>%
    filter(n() > 1) %>% summarise(duplicated_clones = list(clone)) %>% dplyr::pull(duplicated_clones)
  
  #merge if there are duplicated clones
  if (length(duplicated.list) > 0) {
    for (d in 1:length(duplicated.list)){
      rownames(df_uq)[rownames(df_uq) ==  intersect(rownames(df_uq), rownames(df)[seq_along(rownames(df)) %in% duplicated.list[[d]]]) ] <- paste0("s0", d)
      hallmarks$cnv_cluster[hallmarks$cnv_cluster %in% rownames(df)[seq_along(rownames(df)) %in% duplicated.list[[d]]] ] <- paste0("s0", d)
    }
  }
  
  #compute jaccard index and hallmark difference across all pair of clones
  df_var <- as.data.frame(matrix(ncol = 16, nrow = 0))
  labels <- rownames(df_uq)
  for (i in 1:(length(labels)-1)) {
    for (j in (i+1):length(labels)) {
      
      intersect <- sapply(1:ncol(df_uq), function(x) {
        df_uq[labels[i],x] == df_uq[labels[j],x]
      }) %>% sum() +  n_genes[sample]-length(unique(data$gene))
      h_diff <- sapply(paste0("H", 1:13), function(h) {
        mean(hallmarks[hallmarks$cnv_cluster == labels[i], h], na.rm=TRUE) - mean(hallmarks[hallmarks$cnv_cluster == labels[j], h], na.rm=TRUE)
      })
      df_var <- rbind(df_var, c(h_diff, intersect/n_genes[sample], paste0(labels[i], "_", labels[j]), sample))
    }
  }
  df_var[,1:14] <- apply(df_var[,1:14], 2, as.numeric)
  colnames(df_var) <- c(paste0("H", 1:13), "var", "group", "sample")
  
  ## here we will filter/collapse the clones that are very similar
  #first, we identify those clones that are similar and we create a list of similar clones
  low_var <- df_var %>% filter(var >= 0.99)
  low_pair <- list()
  for(pair in low_var$group) {
    p1 <- str_split(pair, "_")[[1]][1]
    p2 <- str_split(pair, "_")[[1]][2]
    if (length(low_pair) == 0) {
      low_pair[[1]] <- c(p1,p2)
    } else {
      found = F
      for (i in 1:length(low_pair)) {
        if (p1 %in% low_pair[[i]] | p2 %in% low_pair[[i]]) {
          low_pair[[i]] <- c(low_pair[[i]], p1, p2)
          found=T
        }
      }
      if (!found) {
        low_pair[[i+1]] <- c(p1,p2)
      }
    }
  }
  low_pair <- lapply(low_pair, unique)
  
  if (length(low_pair) > 1) {
    int <- data.frame(matrix(ncol = 2, nrow = 0))
    for (i in 1:(length(low_pair)-1)) {
      for (j in (i+1):length(low_pair)) {
        if (length(intersect(low_pair[[i]], low_pair[[j]])) > 0) {
          if (!i %in% int[,2]) {
            int <- rbind(int, c(i,j))
          }
        }
      }
    }  
    low_pair_tmp <- low_pair
    for (i in unique(int[,1])) {
      low_pair_tmp[[i]] <- unique(c(low_pair[[i]], do.call(c, lapply(int[int$X1L == i,2], function(x){low_pair[[x]]}))))
    }
    for (i in sort(unique(int[,2]), decreasing = T)) {
      low_pair_tmp[[i]] <- NULL
      print(i)
    }
    low_pair <- low_pair_tmp
  }
  #this function gets the mode of a vector
  get_mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  #then, we will collapse those pair clones that are similar and create the CNV profile by getting the mode
  if (length(low_pair) > 0) {
    df_uq_new <- df_uq[!rownames(df_uq) %in% do.call(c, low_pair),] 
    for (i in 1:length(low_pair)) {
      collapsed <- df_uq[low_pair[[i]],] %>% apply(2,get_mode)
      df_uq_new <- rbind(df_uq_new, collapsed)
      rownames(df_uq_new)[nrow(df_uq_new)] <- paste0("m",i)
      
      hallmarks$cnv_cluster[hallmarks$cnv_cluster %in% low_pair[[i]]  ] <- paste0("m",i)
    }
    df_var <- as.data.frame(matrix(ncol = 16, nrow = 0))
    labels <- rownames(df_uq_new)
    
    #compute jaccard index across all pair of clones after collapsing clones
    for (i in 1:(length(labels)-1)) {
      for (j in (i+1):length(labels)) {
        
        intersect <- sapply(1:ncol(df_uq_new), function(x) {
          df_uq_new[labels[i],x] == df_uq_new[labels[j],x]
        }) %>% sum() +  n_genes[sample]-length(unique(data$gene))
        h_diff <- sapply(paste0("H", 1:13), function(h) {
          mean(hallmarks[hallmarks$cnv_cluster == labels[i], h], na.rm=TRUE) - mean(hallmarks[hallmarks$cnv_cluster == labels[j], h], na.rm=TRUE)
        })
        df_var <- rbind(df_var, c(h_diff, intersect/n_genes[sample], paste0(labels[i], "_", labels[j]), sample))
      }
    }
    df_var[,1:14] <- apply(df_var[,1:14], 2, as.numeric)
    colnames(df_var) <- c(paste0("H", 1:13), "var", "group", "sample")  
    
  }
  #store the hallmarks scores with the filter clones labels
  write.table(hallmarks[,c("cnv_cluster", "cnv_cluster2")] %>% unique(), paste0("/df/", sample, "_map_bayes_filt.txt"))
  #store the pair of clones with genomic distance (jaccard index) and hallmark differences
  write.table(df_var, paste0("/df/", sample, "_jaccard_bayes_filt.txt"))
}

