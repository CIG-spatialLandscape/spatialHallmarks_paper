


samples <- c("Breast", "BreastA", "cHC-1T", "Colorectal", "CRC1",
  "CRC2",  "Ductal", "Ductal_FFPE", "Glioblastoma","ffpe_c_7", "ffpe_c_20", "ffpe_c_21",
  "ffpe_c_34", "ffpe_c_51", "HCC-1T", "HCC-2T", "HCC-5D",
   "ICC-1L", "Intestine", "OV4A", "Ovarian", "Ovarian_FFPE", "OVD1",
  "TNBC-A", "UKF242-T", "UKF260-T", "UKF269-T", "UKF275-T")

#missing Acinar, OV4A, IC, CUP295, PDAC

DU <- c("DU2", "DU3", "DU8", "DU12", "DU13")

"Desktop/IJC/datasets/Public/Breast/RDS/"
#check if files exists

for (sample in samples) {
  if (!file.exists(paste0("Desktop/IJC/datasets/Public/", sample, "/RDS/", sample, ".rds"))) print(sample)
}

for (sample in DU) {
  if (!file.exists(paste0("Desktop/IJC/datasets/Dutreneo/RDS/Assays/", sample, ".rds"))) print(sample)
}

df <- data.frame()

for (sample in samples){
  STobject <- readRDS(paste0("Desktop/IJC/datasets/Public/", sample, "/RDS/", sample, ".rds"))
  print(sample)
  if (!is.null(STobject@assays$Spatial)) df <- rbind(df, c(sample, nrow(STobject), ncol(STobject), median(colSums(STobject@assays$Spatial@counts))))
  if (!is.null(STobject@assays$RNA)) df <- rbind(df, c(sample, nrow(STobject), ncol(STobject), median(colSums(STobject@assays$RNA@counts))))
  
  rm(STobject)
  gc()
}


for (sample in DU){
  STobject <- readRDS(paste0("Desktop/IJC/datasets/Dutreneo/RDS/Assays/", sample, ".rds"))
  df <- rbind(df, c(sample, nrow(STobject), ncol(STobject), median(colSums(STobject@assays$RNA@counts))))
}

colnames(df) <- c("Sample", "Nfeature", "Nspots", "Median UMI count")
