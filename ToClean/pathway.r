####################################
#Wrap up to produce all the figures
####################################

## Libraries
library(Seurat)
library(STutility)
library(ggplot2)
library(gprofiler2)


################################################################################
path_figures <- "Desktop/IJC/datasets/Public/HCV2/figures/"
nfactors <- 5
STobject <- readRDS("Desktop/enhanced/HCV2_enhanced_St_With_Halls.rds")

#### Plot factors

#function to select top loadings by using a threshold (1.645 SD of the mean)
select_top <- function(st, factor) {
  vals <- st@reductions$NMF@feature.loadings[, factor]
  logvals <- log(vals)
  logvals <- logvals[is.finite(logvals)]
  thr <- mean(logvals) + 1.4*sd(logvals)
  names(logvals[logvals > thr])
}


top.genes <- list()
for (factor in 1:nfactors) {
  top.genes[[factor]] <- select_top(STobject, factor)
}

for (source in c("GO:BP", "KEGG", "WP", "HPA", "REAC")) {


#compute gene ontology (GO:BP) for each factor
pathways <- lapply(seq_along(top.genes), function(i) {
  gset <- top.genes[[i]]
  df <- gost(query = gset, organism = "hsapiens", sources = source, evcodes = TRUE, )$result
  #df2 <- gost(query = gset, organism = 'gp__aJd5_EBf2_fFc', evcodes = TRUE)$result     #hallmarks from DBsig
  #df <- rbind(df1, df2)
  if (is.null(df)) return(NULL)
  df$factor <- paste0("factor_", i)
  return(df)
})

pathways <- do.call(rbind, pathways)

pathways$GeneRatio <- pathways$intersection_size/pathways$query_size
#select top 30 (by significance)
pathways.summarized <- pathways %>% 
  group_by(factor) %>%
  top_n(n = 30, wt = -log10(p_value))

#create ggplots for each factor (showing top 30)
p.list <- lapply(1:nfactors, function(i) {
  if (!paste0("factor_", i) %in% pathways.summarized$factor) return(NULL)
  g <- subset(pathways.summarized, factor %in% paste0("factor_", i))
  #g$term_name <- ifelse(nchar(g$term_name) > 40, paste0(substr(g$term_name, start = 1, stop = 40), "..."), g$term_name)
  ggplot() +
    geom_point(data = g, aes(reorder(term_name, -log10(p_value)), -log10(p_value), fill = -log10(p_value), size = GeneRatio), color = "black", shape = 21) +
    coord_flip() +
    facet_grid(~factor) +
    scale_size_continuous(range = c(0.5, 8)) +
    scale_fill_gradientn(colours = viridis::magma(n = 9) %>% rev()) +
    theme_minimal() +
    labs(x = "term", y = "")
})

#load images to the corresponding folder
for (factor in 1:nfactors){
  file_name = paste0(path_figures, "factors", nfactors ,"/enrichment/loadings_", source, "_", factor, ".png")
  ggsave(file_name, plot = p.list[[factor]], height = 7, width = 15, bg = "white")
}

}


#######  Clusters

DefaultAssay(STobject) <- "SCT"
Idents(STobject) <- STobject$NMF_clusters
  
markers <- FindAllMarkers(STobject, only.pos = T, assay = "SCT", slot = "data",
                          logfc.threshold = 0.5, group.by = "NMF_clusters") %>% group_by(cluster)




top_30 <- markers %>% top_n(30, avg_log2FC)

top.genes <- list()
nclusters <- 5
for (cluster in 1:nclusters) {
  top.genes[[as.character(cluster)]] <- top_30[top_30$cluster==cluster,]$gene
}


#compute gene ontology (GO:BP) for each factor

for (source in c("GO:BP", "KEGG", "WP", "HPA", "REAC")) {
pathways <- lapply(seq_along(top.genes), function(i) {
  gset <- top.genes[[i]]
  df <- NULL
  if (length(top.genes[[i]]) !=0) {
    df <- gost(query = gset, organism = "hsapiens", sources = source, evcodes = TRUE, )$result
  }
  #df2 <- gost(query = gset, organism = 'gp__aJd5_EBf2_fFc', evcodes = TRUE)$result     #hallmarks from DBsig
  #df <- rbind(df1, df2)
  if (is.null(df)) return(NULL)
  df$cluster <- paste0("cluster_", i)
  return(df)
})

pathways <- do.call(rbind, pathways)

pathways$GeneRatio <- pathways$intersection_size/pathways$query_size
#select top 30 (by significance)
pathways.summarized <- pathways %>% 
  group_by(cluster) %>%
  top_n(n = 30, wt = -log10(p_value))

#create ggplots for each factor (showing top 30)
p.list <- lapply(1:nclusters, function(i) {
  if (!paste0("cluster_", i) %in% pathways.summarized$cluster) return(NULL)
  g <- subset(pathways.summarized, cluster %in% paste0("cluster_", i))
  #g$term_name <- ifelse(nchar(g$term_name) > 40, paste0(s ubstr(g$term_name, start = 1, stop = 40), "..."), g$term_name)
  ggplot() +
    geom_point(data = g, aes(reorder(term_name, -log10(p_value)), -log10(p_value), fill = -log10(p_value), size = GeneRatio), color = "black", shape = 21) +
    coord_flip() +
    facet_grid(~cluster) +
    scale_size_continuous(range = c(0.5, 8)) +
    scale_fill_gradientn(colours = viridis::magma(n = 9) %>% rev()) +
    theme_minimal() +
    labs(x = "term", y = "")
})

#load images to the corresponding folder
for (cluster in 1:nclusters){
  file_name = paste0(path_figures, "factors5/enrichment/DEG_", source, "_", cluster, ".png")
  ggsave(file_name, plot = p.list[[cluster]], height = 7, width = 15, bg = "white")
}


}







files <- list.files(path="Desktop/enhanced", pattern="*.rds", full.names=TRUE, recursive=FALSE)

for (file in files) {
  name <- strsplit(strsplit(file, split = "/")[[1]][3], split = "_")[[1]][1]
  path_figures <- paste0("Desktop/IJC/datasets/Public/",name, "/figures/") 
  nfactors <- 5
  STobject <- readRDS(file)
  
  #### Plot factors
  
  #function to select top loadings by using a threshold (1.645 SD of the mean)
  
  sdval <- 1.4
  if (name %in% c("Acinar", "Ovarian")) sdval <- 1.1
  
  select_top <- function(st, factor) {
    vals <- st@reductions$NMF@feature.loadings[, factor]
    logvals <- log(vals)
    logvals <- logvals[is.finite(logvals)]
    thr <- mean(logvals) + sdval*sd(logvals)
    names(logvals[logvals > thr])
  }
  
  
  top.genes <- list()
  for (factor in 1:nfactors) {
    top.genes[[factor]] <- select_top(STobject, factor)
  }
  
  for (source in c("GO:BP", "KEGG", "WP", "HPA", "REAC")) {
    
    
    #compute gene ontology (GO:BP) for each factor
    pathways <- lapply(seq_along(top.genes), function(i) {
      gset <- top.genes[[i]]
      df <- gost(query = gset, organism = "hsapiens", sources = source, evcodes = TRUE, )$result
      #df2 <- gost(query = gset, organism = 'gp__aJd5_EBf2_fFc', evcodes = TRUE)$result     #hallmarks from DBsig
      #df <- rbind(df1, df2)
      if (is.null(df)) return(NULL)
      df$factor <- paste0("factor_", i)
      return(df)
    })
    
    pathways <- do.call(rbind, pathways)
    
    pathways$GeneRatio <- pathways$intersection_size/pathways$query_size
    #select top 30 (by significance)
    pathways.summarized <- pathways %>% 
      group_by(factor) %>%
      top_n(n = 30, wt = -log10(p_value))
    
    #create ggplots for each factor (showing top 30)
    p.list <- lapply(1:nfactors, function(i) {
      if (!paste0("factor_", i) %in% pathways.summarized$factor) return(NULL)
      g <- subset(pathways.summarized, factor %in% paste0("factor_", i))
      #g$term_name <- ifelse(nchar(g$term_name) > 40, paste0(substr(g$term_name, start = 1, stop = 40), "..."), g$term_name)
      ggplot() +
        geom_point(data = g, aes(reorder(term_name, -log10(p_value)), -log10(p_value), fill = -log10(p_value), size = GeneRatio), color = "black", shape = 21) +
        coord_flip() +
        facet_grid(~factor) +
        scale_size_continuous(range = c(0.5, 8)) +
        scale_fill_gradientn(colours = viridis::magma(n = 9) %>% rev()) +
        theme_minimal() +
        labs(x = "term", y = "")
    })
    
    #load images to the corresponding folder
    for (factor in 1:nfactors){
      file_name = paste0(path_figures, "factors", nfactors ,"/enrichment/loadings_", source, "_", factor, ".png")
      ggsave(file_name, plot = p.list[[factor]], height = 7, width = 15, bg = "white")
    }
    
  }
  
  
  DefaultAssay(STobject) <- "SCT"
  Idents(STobject) <- STobject$NMF_clusters
  
  markers <- FindAllMarkers(STobject, only.pos = T, assay = "SCT", slot = "data",
                            logfc.threshold = 0.5, group.by = "NMF_clusters") %>% group_by(cluster)
  gc()
  top_30 <- markers %>% top_n(30, avg_log2FC)
  
  top.genes <- list()
  nclusters <- 5
  for (cluster in 1:nclusters) {
    top.genes[[as.character(cluster)]] <- top_30[top_30$cluster==cluster,]$gene
  }
  
  
  #compute gene ontology (GO:BP) for each factor
  
  for (source in c("GO:BP", "KEGG", "WP", "HPA", "REAC")) {
    pathways <- lapply(seq_along(top.genes), function(i) {
      gset <- top.genes[[i]]
      df <- NULL
      if (length(top.genes[[i]]) !=0) {
        df <- gost(query = gset, organism = "hsapiens", sources = source, evcodes = TRUE, )$result
      }
      #df2 <- gost(query = gset, organism = 'gp__aJd5_EBf2_fFc', evcodes = TRUE)$result     #hallmarks from DBsig
      #df <- rbind(df1, df2)
      if (is.null(df)) return(NULL)
      df$cluster <- paste0("cluster_", i)
      return(df)
    })
    
    pathways <- do.call(rbind, pathways)
    
    pathways$GeneRatio <- pathways$intersection_size/pathways$query_size
    #select top 30 (by significance)
    pathways.summarized <- pathways %>% 
      group_by(cluster) %>%
      top_n(n = 30, wt = -log10(p_value))
    
    #create ggplots for each factor (showing top 30)
    p.list <- lapply(1:nclusters, function(i) {
      if (!paste0("cluster_", i) %in% pathways.summarized$cluster) return(NULL)
      g <- subset(pathways.summarized, cluster %in% paste0("cluster_", i))
      #g$term_name <- ifelse(nchar(g$term_name) > 40, paste0(s ubstr(g$term_name, start = 1, stop = 40), "..."), g$term_name)
      ggplot() +
        geom_point(data = g, aes(reorder(term_name, -log10(p_value)), -log10(p_value), fill = -log10(p_value), size = GeneRatio), color = "black", shape = 21) +
        coord_flip() +
        facet_grid(~cluster) +
        scale_size_continuous(range = c(0.5, 8)) +
        scale_fill_gradientn(colours = viridis::magma(n = 9) %>% rev()) +
        theme_minimal() +
        labs(x = "term", y = "")
    })
    
    #load images to the corresponding folder
    for (cluster in 1:nclusters){
      file_name = paste0(path_figures, "factors5/enrichment/DEG_", source, "_", cluster, ".png")
      ggsave(file_name, plot = p.list[[cluster]], height = 7, width = 15, bg = "white")
    }
    
    
  }
  rm(STobject)
  gc()
}




