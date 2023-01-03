

library(Seurat)
library(dplyr)

supp_tbl <- data.frame(matrix(ncol=10, nrow=1))
colnames(supp_tbl) <- c("Tumor type", "Hallmark", "Factor","Number of top contributing genes", 
                        "Number of hallmark genes found in normalized ST data",
                        "Percentage of hallmarks genes in loading genes",
                        "Percentage of loading genes in hallmarks genes", 
                        "Number of top marker genes (DEA)", 
                        "Percentage of hallmarks genes in DEA genes",
                        "Percentage of DEA genes in hallmarks genes" 
                        )
files <- list.files(path="Desktop/enhanced", pattern="*.rds", full.names=TRUE, recursive=FALSE)
files <- files[-c(6,7,8,10)]
files <- files[-c(1,2,3)]
for (file in files) {
  STobject <- readRDS(file)
  name <- strsplit(strsplit(file, "/")[[1]][3], "_")[[1]][1] 
  print(name)
  
  sd_val <- 1.4
  if (name %in% c("Ovarian", "IC", "Acinar")) sd_val <- 1.1
  select_top <- function(st, factor) {
    vals <- st@reductions$NMF@feature.loadings[, factor]
    logvals <- log(vals)
    logvals <- logvals[is.finite(logvals)]
    thr <- mean(logvals) + sd_val*sd(logvals)
    names(logvals[logvals > thr])
  }
  
  
  
  DefaultAssay(STobject) <- "SCT"
  Idents(STobject) <- STobject$NMF_clusters
  markers <- FindAllMarkers(STobject, only.pos = T, assay = "SCT", slot = "data",
                            logfc.threshold = 0.5, group.by = "NMF_clusters") %>% group_by(cluster)
  top_30 <- markers %>% top_n(30, avg_log2FC)
  gc()
  
  colnames(top_30)[6] <- "factor" 
  top_30$`Tumor type` <- name
  write.table(top_30, paste0("Desktop/DEA/", name, ".txt"), sep = "\t", row.names = F)
  
  for (factor in 1:5) {
    loadings <- select_top(STobject, factor)
     DEA <- top_30[top_30$factor==factor,]$gene
      for(hallmark in 1:13) {
        hallmark_genes <- H_features[[hallmark]][H_features[[hallmark]] %in% rownames(STobject[["SCT"]])]
        tbl_loadings <- table(hallmark_genes %in% loadings)
        tbl2_loadings <- table(loadings %in% hallmark_genes)
        
        
        tbl_DEA <- table(hallmark_genes %in% DEA)
        tbl2_DEA <- table(DEA %in% hallmark_genes)
        
        
        supp_tbl <- rbind(supp_tbl, c(name, hallmark_names[hallmark], factor, 
                                      length(loadings), length(hallmark_genes), as.numeric(tbl_loadings["TRUE"])/as.numeric(sum(tbl_loadings))*100,
                                      as.numeric(tbl2_loadings["TRUE"])/as.numeric(sum(tbl2_loadings))*100,
                                      length(DEA), as.numeric(tbl_DEA["TRUE"])/as.numeric(sum(tbl_DEA))*100,
                                      as.numeric(tbl2_DEA["TRUE"])/as.numeric(sum(tbl2_DEA))*100))
      }
  }
  rm(STobject)
  gc()
}

supp_tbl <- supp_tbl[-1,]

supp_tbl[is.na(supp_tbl)] <- 0
supp_tbl$`Percentage of loading genes in hallmarks genes` <- as.numeric(supp_tbl$`Percentage of loading genes in hallmarks genes`)
supp_tbl$`Percentage of hallmarks genes in loading genes` <- as.numeric(supp_tbl$`Percentage of hallmarks genes in loading genes`)
supp_tbl$`Percentage of DEA genes in hallmarks genes` <- as.numeric(supp_tbl$`Percentage of DEA genes in hallmarks genes`)
supp_tbl$`Percentage of hallmarks genes in DEA genes` <- as.numeric(supp_tbl$`Percentage of hallmarks genes in DEA genes`)


supp_tbl$`Tumor type`[supp_tbl$`Tumor type`=="Breast"] <- "Invasive Lobular Carcinoma (Breast)"
supp_tbl$`Tumor type`[supp_tbl$`Tumor type`=="Ductal"] <- "Invasive Ductal Carcinoma (Breast)"
supp_tbl$`Tumor type`[supp_tbl$`Tumor type`=="IC"] <- "Adenocarcinoma (Prostate)"
supp_tbl$`Tumor type`[supp_tbl$`Tumor type`=="OV4A"] <- "HGSOC (Ovary)"
supp_tbl$`Tumor type`[supp_tbl$`Tumor type`=="Ovarian"] <- "Endometrial Adenocarcinoma (Ovary)"
supp_tbl$`Tumor type`[supp_tbl$`Tumor type`=="Acinar"] <- "Acinar Cell Carcinoma (Prostate)"
supp_tbl$`Tumor type`[supp_tbl$`Tumor type`=="Colorectal"] <- "Invasive Adenocarcinoma (Colorectal)"

write.table(supp_tbl, "Desktop/IJC/datasets/IGTP/figuresPaper/tables/tbl_loadings.xls", sep = "\t", row.names = F)


x <- read.table("Desktop/IJC/datasets/IGTP/figuresPaper/tables/tbl_S2.xls", header = T, sep = "\t")
x$Tumor.type[x$Tumor.type=="OV4A"] <- "HGSOC (Ovary)"

x <- write.table(x, "Desktop/IJC/datasets/IGTP/figuresPaper/tables/tbl_S2.xls", row.names = F, sep = "\t")

supp_tbl2 <- data.frame(matrix(ncol=9, nrow=1))
colnames(supp_tbl2) <- c("Tumor type", "Hallmark", "Factor","Number of intersection genes (loadings)",
                        "Percentage of loading genes in hallmark genes (intersection)", "Number of intersection genes (DEA)",
                        "Percentage of DEA genes in hallmark genes (intersection)", "Percentage of hallmarks genes in loading genes (intersection)",
                        "Percentage of hallmarks genes in DEA genes (intersection)")
for (tumor in unique(x$Tumor.type)) {
  for (factor in unique(x$factor)) {
    genes_loadings <- unique(do.call(c, strsplit(x[x$Tumor.type==tumor & x$factor==factor & x$method=="loadings",]$intersection, split = ",")))
    genes_DEA <- unique(do.call(c, strsplit(x[x$Tumor.type==tumor & x$factor==factor & x$method=="DEA",]$intersection, split = ",")))
     for(hallmark in 1:13) {
       tbl_loadings <- table(genes_loadings %in% H_features[[hallmark]])
       tbl_DEA <- table(genes_DEA %in% H_features[[hallmark]])
       supp_tbl2 <- rbind(supp_tbl2, c(tumor, hallmark_names[hallmark], factor, 
                                  length(genes_loadings), as.numeric(tbl_loadings["TRUE"])/as.numeric(sum(tbl_loadings))*100,
                                  length(genes_DEA), as.numeric(tbl_DEA["TRUE"])/as.numeric(sum(tbl_DEA))*100, 
                          as.numeric(tbl_loadings["TRUE"]), as.numeric(tbl_DEA["TRUE"])))
     }
   }
}
supp_tbl2 <- supp_tbl2[-1,]
supp_tbl2[is.na(supp_tbl2)] <- 0

write.table(supp_tbl2, "Desktop/IJC/datasets/IGTP/figuresPaper/tables/tbl_intersection.xls", sep = "\t", row.names = F)

supp_tbl$Factor <- paste0("factor_", supp_tbl$Factor)













supp_tbl <- data.frame(matrix(ncol=10, nrow=1))
colnames(supp_tbl) <- c("Tumor type", "Hallmark", "Factor","Number of top contributing genes", 
                        "Number of hallmark genes found in normalized ST data",
                        "Percentage of loading genes in hallmarks genes", 
                        "Percentage of hallmarks genes in loading genes", 
                        "Number of top marker genes (DEA)", 
                        "Percentage of DEA genes in hallmarks genes", 
                        "Percentage of hallmarks genes in DEA genes")
files <- list.files(path="Desktop/enhanced", pattern="*.rds", full.names=TRUE, recursive=FALSE)
files <- files[-c(6,7,8,10)]
for (file in files) {
  STobject <- readRDS(file)
  name <- strsplit(strsplit(file, "/")[[1]][3], "_")[[1]][1] 
  print(name)
  
  sd_val <- 1.4
  if (name %in% c("Ovarian", "IC", "Acinar")) sd_val <- 1.1
  select_top <- function(st, factor) {
    vals <- st@reductions$NMF@feature.loadings[, factor]
    logvals <- log(vals)
    logvals <- logvals[is.finite(logvals)]
    thr <- mean(logvals) + sd_val*sd(logvals)
    names(logvals[logvals > thr])
  }
  top_30 <- read.table(paste0("Desktop/DEA/", name, ".txt"), header = T, sep = "\t")
  for (factor in 1:5) {
    loadings <- select_top(STobject, factor)
    DEA <- top_30[top_30$factor==factor,]$gene
    for(hallmark in 1:13) {
      hallmark_genes <- H_features[[hallmark]][H_features[[hallmark]] %in% rownames(STobject[["SCT"]])]
      tbl_loadings <- table(hallmark_genes %in% loadings)
      tbl2_loadings <- table(loadings %in% hallmark_genes)
      
      
      tbl_DEA <- table(hallmark_genes %in% DEA)
      tbl2_DEA <- table(DEA %in% hallmark_genes)
      
      
      supp_tbl <- rbind(supp_tbl, c(name, hallmark_names[hallmark], factor, 
                                    length(loadings), length(hallmark_genes), as.numeric(tbl_loadings["TRUE"])/as.numeric(sum(tbl_loadings))*100,
                                    as.numeric(tbl2_loadings["TRUE"])/as.numeric(sum(tbl2_loadings))*100,
                                    length(DEA), as.numeric(tbl_DEA["TRUE"])/as.numeric(sum(tbl_DEA))*100,
                                    as.numeric(tbl2_DEA["TRUE"])/as.numeric(sum(tbl2_DEA))*100))
    }
  }
  rm(STobject)
  gc()
}

supp_tbl <- supp_tbl[-1,]

supp_tbl[is.na(supp_tbl)] <- 0

write.table(supp_tbl, "Desktop/IJC/datasets/IGTP/figuresPaper/tables/tbl_DEA&loadings.xls", sep = "\t", row.names = F)


tbl <- cbind(supp_tbl, supp_tbl2[,4:7])
colnames(tbl)[c(12,14)] <- c("Percentage of loading genes in hallmarks genes (intersection)", "Percentage of DEA genes in hallmarks genes (intersection)")
tbl[,5] <- as.numeric(tbl[,5])
tbl[,11] <- as.numeric(tbl[,11])
tbl[,12] <- as.numeric(tbl[,12])
tbl[,13] <- as.numeric(tbl[,13])
tbl[,14] <- as.numeric(tbl[,14])


supp_tbl2$`Percentage of hallmarks genes in loading genes (intersection)` <- as.numeric(tbl[,8])/tbl[,5]*100
supp_tbl2$`Percentage of hallmarks genes in DEA genes (intersection)` <-  as.numeric(tbl[,9])/tbl[,5]*100

supp_tbl2 <- supp_tbl2[,c(1,2,3,10,4,8,5,6,9,7)]
colnames(supp_tbl2)[4] <- "Number of hallmark genes found in normalized ST data"
supp_tbl2[,4:10] <- apply(supp_tbl2[,4:10], 2, as.numeric)
supp_tbl[,4:10] <- apply(supp_tbl[,4:10], 2, as.numeric)
write.table(supp_tbl2, "Desktop/IJC/datasets/IGTP/figuresPaper/tables/tbl_intersection.xls", sep = "\t", row.names = F)
