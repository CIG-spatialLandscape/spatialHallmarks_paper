
samples <- list.files("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/output_df/") %>% str_replace(".txt", "")


files <- list.files("Desktop/IJC/datasets/IGTP/figuresPaper/intersect/")
df <- c()
for (file in files) {
  tbl <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/intersect/", file))
  df <- rbind(df, c(file %>% str_replace("_genes.txt", ""),  nrow(tbl)))
}


summary <- rbind(c("Acinar", "Prostate", "FFPE", "Public", 3038, 0),
      c("Breast", "Breast", "OCT", "Public", 4318, 0),
      c("BreastA", "Breast", "OCT", "Public", 3798, 16006),
      c("Colorectal", "Colorectal", "OCT", "Public", 3115, 0),
      c("Ductal", "Breast", "OCT", "Public", 4704, 0),
      c("DuctalFFPE", "Breast", "FFPE", "Public", 4315, 0),
      c("Glioblastoma", "Glioblastoma", "OCT", "Public", 3038, 0),
      c("HCC1T", "Liver", "OCT", "Public", 3148, 15243),
      c("IC", "Prostate", "FFPE", "Public", 4363, 0),
      c("ICC1L", "Liver", "OCT", "Public", 3756, 0),
      c("Intestine", "Colorectal", "FFPE", "Public", 2657, 0),
      c("OV4A", "Ovary", "OCT", "Self-Generated", 2619, 0),
      c("Ovarian", "Ovary", "OCT", "Public", 3473, 0),
      c("OVFFPE", "Ovary", "FFPE", "Public", 3448, 0),
      c("cHC1T", "Liver", "OCT", "Public", 4759, 0),
      c("TNBCA", "Breast", "OCT", "Public", 4784, 0),
      c("UKF242T", "Glioblastoma", "OCT", "Public", 3029, 0),
      c("UKF260T", "Glioblastoma", "OCT", "Public", 2997, 0),
      c("UKF269T", "Glioblastoma", "OCT", "Public", 3217, 0),
      c("UKF275T", "Glioblastoma", "OCT", "Public", 3726, 0),
      c("CUP295", "Ovary", "FFPE", "Self-Generated", 3804, 0),
      c("CRC1", "Colorectal", "OCT", "Public", 3299, 0),
      c("CRC2", "Colorectal", "OCT", "Public", 4167, 0),
      c("HCC5D", "Liver", "OCT", "Public", 4350, 14583),
      c("C7", "Kidney", "FFPE", "Public", 4974, 0),
      c("C20", "Kidney", "FFPE", "Public", 4947, 0),
      c("C21", "Kidney", "FFPE", "Public", 4914, 0),
      c("C34", "Kidney", "FFPE", "Public", 3584, 0),
      c("C51", "Kidney", "FFPE", "Public", 4359, 0),
      c("OVD1", "Ovary", "FFPE", "Public", 3663, 14814),
      c("PC1", "Prostate", "OCT", "Public", 2618, 0),
      c("PC2", "Prostate", "OCT", "Public", 1933, 0),
      c("Co1", "Colorectal", "FFPE", "Self-Generated", 4760, 16703),
      c("Co2", "Colorectal", "FFPE", "Self-Generated", 4759, 16657),
      c("Co3", "Colorectal", "FFPE", "Self-Generated", 4728, 16061),
      c("Co4", "Colorectal", "FFPE", "Self-Generated", 3587, 15690),
      c("M1", "Breast", "FFPE", "Self-Generated", 3954, 0),
      c("M2", "Breast", "FFPE", "Self-Generated", 4189, 0),
      c("M3", "Breast", "FFPE", "Self-Generated", 4512, 0),
      c("M4", "Breast", "FFPE", "Self-Generated", 4769, 0),
      c("P1", "Lung", "FFPE", "Self-Generated", 4769, 0),
      c("P3", "Lung", "FFPE", "Self-Generated", 4619, 0),
      c("P4", "Lung", "FFPE", "Self-Generated", 4767, 0),
      c("P5", "Lung", "FFPE", "Self-Generated", 4682, 0),
      c("P6", "Lung", "FFPE", "Self-Generated", 4082, 0),
      c("P7", "Lung", "FFPE", "Self-Generated", 4724, 15452),
      c("P8", "Lung", "FFPE", "Self-Generated", 4562, 0),
      c("DU2", "Bladder", "FFPE", "Self-Generated", 3765, 15692),
      c("DU3", "Bladder", "FFPE", "Self-Generated", 2930, 15642),
      c("DU8", "Bladder", "FFPE", "Self-Generated", 4107, 17716),
      c("DU12", "Bladder", "FFPE", "Self-Generated", 2589, 16311),
      c("DU13", "Bladder", "FFPE", "Self-Generated", 3536, 16642),
      c("HCC2T", "Liver", "OCT", "Public", 4724, 0),
      c("P259_H2A2", "Pancreas", "OCT", "Public", 3927, 0),
      c("P264", "Pancreas", "OCT", "Public", 4447, 0),
      c("P270", "Pancreas", "OCT", "Public", 4770, 0),
      c("P288", "Pancreas", "OCT", "Public", 3928, 0),
      c("P306", "Pancreas", "OCT", "Public", 4779, 0)
)


rownames(summary) <- summary[,1]
summary[df[,1],6] <- df[,2]

samples <- summary[summary[,6]=="0",1]
for (sample in samples) {
  sce <- readRDS(paste0("Desktop/IJC/datasets/Public/", sample, "/RDS/", sample, "_sce.rds"))
  summary[sample, 6] <- nrow(sce)
}



samples <- rownames(summary)
for (sample in samples) {
  if (!file.exists(paste0("Desktop/IJC/datasets/Public/", sample, "/RDS/", sample, ".rds"))){ print(sample)}
}

summary <- cbind(summary, rep(0,58))
for (sample in samples[-45]) {
  print(sample)
  STobject <- readRDS(paste0("Desktop/IJC/datasets/Public/", sample, "/RDS/", sample, ".rds"))
  summary[sample, 7] <- median(apply(as.matrix(STobject@assays[[1]]@counts), 2, function(x) {sum(x>0)}))
  rm(STobject)
  gc()
}
STobject <- readRDS(paste0("Desktop/IJC/datasets/Public/", sample, "/RDS/", sample, ".rds"))
median(apply(as.matrix(STobject@assays$RNA@counts), 2, function(x) {sum(x>0)}))

colnames(summary) <- c("Sample", "Tissue location", "Embedding Method", "Source",
                       "Number of spots", "Number of detected genes", "Median unique genes per spot")
summary <- data.frame(summary)
summary[,5:7] <- apply(summary[,5:7],2, as.numeric)

summary <- summary %>% arrange(Tissue.location, Embedding.Method, Source, -Number.of.spots, -Number.of.detected.genes)

library(ComplexHeatmap)
col_tissue_map <- setNames(Seurat::DiscretePalette(15, palette = "polychrome") [c(1:3, 5, 6:11)], sort(unique(summary$Tissue.location)))

col_ha = HeatmapAnnotation(`Tissue` = summary$Tissue.location,
                           `Embedding Method` = summary$Embedding.Method,
                           `Source` = summary$Source, 
                           col = list(Source=c("Public"="gold2", "Self-Generated"="magenta1"),
                                      Tissue=col_tissue_map,
                                      `Embedding Method`=c("FFPE"=Seurat::DiscretePalette(15, palette = "polychrome")[11],
                                                           "OCT"=Seurat::DiscretePalette(15, palette = "polychrome")[14]))
                           )
col_fun1 = colorRamp2(c(min(summary[, 5]), max(summary[, 5])), c("white", "#E69F00"))
col_fun2 = colorRamp2(c(min(summary[, 6]), max(summary[, 6])), c("white", "#009E73"))
col_fun3 = colorRamp2(c(min(summary[, 7]), max(summary[, 7])), c("white", "#CC79A7"))
h1 <- Heatmap(t(summary[, 5]), cluster_columns = F, 
              top_annotation = col_ha, height =  unit(0.5, "cm"), 
              row_labels = "Number of spots", col = col_fun1, name ="Number of spots",
              heatmap_legend_param = list(nrow=1))
h2 <- Heatmap(t(summary[, 6]), cluster_columns = F, name ="Number of detected genes",
              height =  unit(0.5, "cm"), row_labels = "Number of detected genes", 
              col=col_fun2,
              heatmap_legend_param = list(nrow=1))
h3 <- Heatmap(t(summary[, 7]), cluster_columns = F, name = "Median unique genes per spot",
              height =  unit(0.5, "cm"), row_labels = "Median unique genes per spot", 
              col=col_fun3,heatmap_legend_param = list(
                title = "Median unique genes per spot", at = c(1100, 3000,6000,8200),
                nrow=1
              ))
ht_list = h1 %v% h2 %v% h3
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/Final/supp_summary.pdf", height = 6, width = 10)
draw(ht_list,  ht_gap = unit(0.1, "cm"), legend_grouping = "original",
     heatmap_legend_side = "bottom", merge_legend=F,annotation_legend_side = "bottom" )
dev.off()
