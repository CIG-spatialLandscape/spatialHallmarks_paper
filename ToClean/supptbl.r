STobject <- readRDS("Desktop/enhanced/OV4A_enhanced.rds")




files <- files <- list.files(path="Desktop/enhanced", pattern="*.rds", full.names=TRUE, recursive=FALSE)
files <- files[c(1,2,3,4,5,9,11,12)]
for (file in files) {
  STobject <- readRDS(file)
  name <- strsplit(strsplit(file, "/")[[1]][3], "_")[[1]][1]
  print(name)
  #loadings
  sd_val <- 1.4
  if (name %in% c("IC", "Acinar", "Ovarian")) sd_val <- 1.1
  
  select_top <- function(st, factor) {
    vals <- st@reductions$NMF@feature.loadings[, factor]
    logvals <- log(vals)
    logvals <- logvals[is.finite(logvals)]
    thr <- mean(logvals) + sd_val*sd(logvals)
    names(logvals[logvals > thr])
  }
  
  
  top.genes <- list()
  for (factor in 1:5) {
    top.genes[[factor]] <- select_top(STobject, factor)
  }
  
  all_pathways <- c()
  for (source in c("GO:BP", "KEGG", "WP", "HPA", "REAC")) {
    #compute gene ontology (GO:BP) for each factor
    pathways <- lapply(seq_along(top.genes), function(i) {
      gset <- top.genes[[i]]
      df <- gost(query = gset, organism = "hsapiens", sources = source, evcodes = TRUE, )$result
      #df2 <- gost(query = gset, organism = 'gp__aJd5_EBf2_fFc', evcodes = TRUE)$result     #hallmarks from DBsig
      #df <- rbind(df1, df2)
      if (is.null(df)) return(NULL)
      df$factor <- paste0("factor_", i)
      df$method <- "loadings"
      return(df)
    })
    
    pathways <- do.call(rbind, pathways)
    
    pathways$GeneRatio <- pathways$intersection_size/pathways$query_size
    
    #select top 30 (by significance)
    pathways.summarized <- pathways %>% 
      group_by(factor) %>%
      top_n(n = 20, wt = -log10(p_value))
    all_pathways <- rbind(all_pathways, pathways.summarized)
  }
  
  #DEA
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
      df$factor <- paste0("factor_", i)
      df$method <- "DEA"
      return(df)
    })
    
    pathways <- do.call(rbind, pathways)
    
    pathways$GeneRatio <- pathways$intersection_size/pathways$query_size
    #select top 30 (by significance)
    pathways.summarized <- pathways %>% 
      group_by(factor) %>%
      top_n(n = 20, wt = -log10(p_value))
    
    all_pathways <- rbind(all_pathways, pathways.summarized)
  }
  
  all_pathways$CancerType <- name
  all_pathways <- apply(as.data.frame(all_pathways),2,as.character)
  write.table(all_pathways, paste0("Desktop/IJC/datasets/IGTP/figuresPaper/supp_tbl/", name, ".txt"), sep = "\t")
  rm(STobject)
  gc()
}

files <- list.files(path="Desktop/IJC/datasets/IGTP/figuresPaper/supp_tbl", pattern="*.txt", full.names=TRUE, recursive=FALSE)
files <- files[-2]
list.tbl <- lapply(files, function(x){ 
  read.table(x, sep = "\t", header = T)})

all <- do.call(rbind, list.tbl)



all$Cancer[all$Cancer=="Breast"] <- "Invasive Lobular Carcinoma (Breast)"
all$Cancer[all$Cancer=="Ductal"] <- "Invasive Ductal Carcinoma (Breast)"
all$Cancer[all$Cancer=="IC"] <- "Adenocarcinoma (Prostate)"
all$Cancer[all$Cancer=="OV_4A"] <- "HGSOC (Ovary)"
all$Cancer[all$Cancer=="Ovarian"] <- "Endometrial Adenocarcinoma (Ovary)"
all$Cancer[all$Cancer=="Acinar"] <- "Acinar Cell Carcinoma (Prostate)"
all$Cancer[all$Cancer=="Colorectal"] <- "Invasive Adenocarcinoma (Colorectal)"

write.table(all, "Desktop/IJC/datasets/IGTP/figuresPaper/supp_tbl/all.xls", sep = "\t")




names <- c("Acinar", "IC", "Breast", "Ductal", "OV4A", "Ovarian", "Colorectal", "Glioblastoma")


list.stats <- lapply(names, function(x){ 
  stats <- read.table(paste0("Desktop/IJC/datasets/Public/",x,"/figures/hallmarks/stats.txt"), sep = "\t", header = T)
  stats$hallmark <- rownames(stats)
  stats$Cancer <- name
  stats
})
all_stats <- do.call(rbind, list.stats)
rownames(all_stats) <- NULL


all_stats$Cancer[all_stats$Cancer=="Breast"] <- "Invasive Lobular Carcinoma (Breast)"
all_stats$Cancer[all_stats$Cancer=="Ductal"] <- "Invasive Ductal Carcinoma (Breast)"
all_stats$Cancer[all_stats$Cancer=="IC"] <- "Adenocarcinoma (Prostate)"
all_stats$Cancer[all_stats$Cancer=="OV_4A"] <- "HGSOC (Ovary)"
all_stats$Cancer[all_stats$Cancer=="Ovarian"] <- "Endometrial Adenocarcinoma (Ovary)"
all_stats$Cancer[all_stats$Cancer=="Acinar"] <- "Acinar Cell Carcinoma (Prostate)"
all_stats$Cancer[all_stats$Cancer=="Colorectal"] <- "Invasive Adenocarcinoma (Colorectal)"

write.table(all, "Desktop/IJC/datasets/IGTP/figuresPaper/hallmarks_stats.xls", sep = "\t")
