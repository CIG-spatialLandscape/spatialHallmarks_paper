


files <- list.files("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/output_df", full.names = F)
files <- stringr::str_remove(files, pattern = ".txt")


correlations <- list()
for (sample in files) {
  hallmarks <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/objects_mts/", sample, "_hallmarks.txt"), sep = "\t", header = T)
  hallmarks <- hallmarks[rownames(hallmarks)!="subspot_1.1",]
  correlations[[sample]] <- cor(hallmarks)
}

avg <- as.data.frame(do.call(rbind, correlations))
avg$h_names <- rep(paste0("H", 1:13), 58)
avg <- avg %>% group_by(h_names) %>% summarise(H1 = mean(H1),
                                        H2 = mean(H2),
                                        H3 = mean(H3),
                                        H4 = mean(H4),
                                        H5 = mean(H5),
                                        H6 = mean(H6),
                                        H7 = mean(H7),
                                        H8 = mean(H8),
                                        H9 = mean(H9),
                                        H10 = mean(H10),
                                        H11 = mean(H11),
                                        H12 = mean(H12),
                                        H13 = mean(H13)) %>% as.data.frame()
rownames(avg) <- avg$h_names
avg$h_names <- NULL

library(pheatmap)
png(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/correlation_plots/Pancancer_cor.png"), width = 350, height = 350, )
pheatmap(avg, scale = "none", )
dev.off()

###

tumor_type <- list()
tumor_type[["Prostate"]] <- c("Acinar", "IC", "PC1", "PC2")
tumor_type[["Breast"]] <- c("Breast", "BreastA", "Ductal", "DuctalFFPE", "TNBCA", "M1", "M2", "M3", "M4")
tumor_type[["Kidney"]]  <- c("C20", "C21", "C34", "C51", "C7")
tumor_type[["Liver"]] <- c("cHC1T", "HCC1T", "HCC2T", "HCC5D", "ICC1L")
tumor_type[["Colorectal"]] <- c("Colorectal", "CRC1", "CRC2", "Intestine", "Co1", "Co2", "Co3", "Co4")
tumor_type[["Bladder"]] <- c("DU2", "DU3", "DU12", "DU8", "DU13")
tumor_type[["Ovarian"]] <- c("OV4A", "OVD1", "Ovarian", "CUP295", "OVFFPE")
tumor_type[["Pancreas"]] <- c("P259_H2A2", "P264", "P270", "P288", "P306")
tumor_type[["Glioblastoma"]] <- c("Glioblastoma", "UKF242T", "UKF260T", "UKF269T", "UKF275T")
tumor_type[["Lung"]] <- c("P1","P3", "P4", "P5", "P6", "P7","P8")


# Prostate
for (tumor in names(tumor_type)) {
  correlations <- list()
  for (sample in tumor_type[[tumor]]) {
    hallmarks <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/objects_mts/", sample, "_hallmarks.txt"), sep = "\t", header = T)
    hallmarks <- hallmarks[rownames(hallmarks)!="subspot_1.1",]
    correlations[[sample]] <- cor(hallmarks)
  }
  avg <- as.data.frame(do.call(rbind, correlations))
  avg$h_names <- rep(paste0("H", 1:13), length(tumor_type[[tumor]]))
  avg <- avg %>% group_by(h_names) %>% summarise(H1 = mean(H1),
                                                 H2 = mean(H2),
                                                 H3 = mean(H3),
                                                 H4 = mean(H4),
                                                 H5 = mean(H5),
                                                 H6 = mean(H6),
                                                 H7 = mean(H7),
                                                 H8 = mean(H8),
                                                 H9 = mean(H9),
                                                 H10 = mean(H10),
                                                 H11 = mean(H11),
                                                 H12 = mean(H12),
                                                 H13 = mean(H13)) %>% as.data.frame()
  rownames(avg) <- avg$h_names
  avg$h_names <- NULL
  png(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/correlation_plots/", tumor, "_cor.png"), width = 350, height = 350, )
  pheatmap(avg, scale = "none", )
  dev.off()
}


