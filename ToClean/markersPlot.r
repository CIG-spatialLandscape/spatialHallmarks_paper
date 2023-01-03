source("Desktop/IJC/datasets/IGTP/figuresPaper/scripts/utilities/BayesSpace_functions.r")
source("Desktop/IJC/datasets/IGTP/figuresPaper/scripts/utilities/AnnotateCorrelations.r")

files <- list.files("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/output_df", full.names = T)
df_all <- lapply(files, function(x) {
  read.table(x, sep = "\t")
})
df_all <- do.call(rbind, df_all)

pdac <- c("P259_H2A2", "P264", "P270", "P288", "P306")
bladder <- c("DU2", "DU3", "DU8", "DU12", "DU13")
colorectal <- c("Colorectal", "Intestine", "CRC1", "CRC2")
liver <- c("cHC1T", "HCC1T", "HCC2T", "HCC5D", "ICC1L")
glioblastoma <- c("Glioblastoma", "UKF242T", "UKF260T", "UKF269T", "UKF275T")
breast <- c("Breast", "BreastA", "Ductal", "DuctalFFPE", "TNBCA")
ovarian <- c("OV4A", "Ovarian", "OVD1", "OVFFPE")
kidney <- c("C7", "C20", "C21", "C34", "C51")
prostate <- c("Acinar", "IC")

samples <- unique(sapply(list.files("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/objects_mts/"), function(x){
  strsplit(x, split = "_")[[1]][1]
}))
samples <- c(samples[1:31], samples[37:41])
#samples[32] <- "P259_H2A2"


#plot Hallmarks
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/enhanced_sce//", sample, "_sce_enhanced.rds"))
  sce <- sce[,-1]
  for(h in paste0("H", 1:13)) {
    colData(sce)[,h] <- df_all[df_all$sample==sample, h]
    
    v <- .make_triangle_subspots(colData(sce), fill = h)
    p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
      scale_fill_viridis_c(option = "B") + labs(fill = "")
    pdf(paste0("Desktop/markers/Hallmarks/", h, "/", sample, "_", h, ".pdf"))
    plot(p)
    dev.off()
  }
  rm(sce)
  gc()
}

#plot Hallmarks
for (sample in pdac) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  for(h in paste0("H", 1:13)) {
    colData(sce)[,h] <- df_all[df_all$sample==sample, h]
    
    v <- .make_triangle_subspots(colData(sce), fill = h)
    p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
      scale_fill_viridis_c(option = "B") + labs(fill = "")
    pdf(paste0("Desktop/markers/Hallmarks/", h, "/", sample, "_", h, ".pdf"))
    plot(p)
    dev.off()
  }
  rm(sce)
  gc()
}


### H1

#EGFR: Pancreas!!!!!!!!!!!!!!, Prostate
for (sample in prostate) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)["EGFR",]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H1/EGFR/",  sample, "_EGFR.pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}

#AKT1: Kidney
gene_name <- "AKT1"
for (sample in kidney) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H1/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}

### H2
#CDK1: Colorectal, Breast, Ovary, Liver, Bladder, Pancreas, Kidney
gene_name <- "CDK1"
samples <- c(breast, ovarian, liver, bladder, kidney)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H2/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}


### H5
#NFKBIA: Prostate, Colorectal
gene_name <- "NFKBIA"
samples <- c(colorectal, prostate)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H5/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}

### H6
#VCAN: Colorectal, Pancreas, Ovary, Bladder, Breast, Liver, kidney
gene_name <- "VCAN"
samples <- c(colorectal, ovarian, bladder, breast, liver, kidney)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H6/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}

### H8
#BRCA1: Colorectal
gene_name <- "BRCA1"
samples <- c(colorectal)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H8/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}
#XRCC4: Liver, pancreas
gene_name <- "XRCC4"
samples <- c(liver)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H8/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}
#TP53BP1: Liver, prostate
gene_name <- "TP53BP1"
samples <- c(liver, prostate)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H8/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}
#MSH2: colorectal, liver, ovarian, breast, prostate, kidney
gene_name <- "MSH2"
samples <- c(colorectal, liver, ovarian, breast, prostate, kidney)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H8/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}
#PMS2: breast, prostate, liver, kidney, colorectal
gene_name <- "PMS2"
samples <- c(colorectal, liver, ovarian, breast, prostate, kidney)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H8/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}

### H10
#ATP6V1B2: brain
gene_name <- "ATP6V1B2"
samples <- c(glioblastoma)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H10/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}
#ATP6V1H: brain
gene_name <- "ATP6V1H"
samples <- c(glioblastoma)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H10/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}
#ATP6V1D: brain
gene_name <- "ATP6V1D"
samples <- c(glioblastoma)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H10/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}
#ATP6V1A: brain
gene_name <- "ATP6V1A"
samples <- c(glioblastoma)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H10/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}
#GAPDH: brain, pancreas
gene_name <- "GAPDH"
samples <- c(glioblastoma)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H10/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}
#VDAC1: brain
gene_name <- "VDAC1"
samples <- c(glioblastoma)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H10/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}
### H11
#TXN: pancreas, colorectal, breast
gene_name <- "TXN"
samples <- c(colorectal, breast)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H11/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}
### H12
#EZH2: liver, pancreas, prostate, breast, kidney, colorectal
gene_name <- "EZH2"
samples <- c(liver, prostate, breast, kidney, colorectal)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H12/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}
#SUZ12: liver, pancreas, prostate, breast, kidney, colorectal
gene_name <- "SUZ12"
samples <- c(liver, prostate, breast, kidney, colorectal)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H12/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}
### H3
#LAG3: kidney
gene_name <- "LAG3"
samples <- c(kidney)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H3/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}
### H4
#TERT: Liver
gene_name <- "TERT"
samples <- c(liver)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H4/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}
#TRF1: Liver
gene_name <- "TRF1"
samples <- c(liver)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H4/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}
#TRF2: Liver
gene_name <- "TRF2"
samples <- c(liver)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H4/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}
#POT1: Liver
gene_name <- "POT1"
samples <- c(liver, prostate)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H4/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}
### H7
#NOTCH3
gene_name <- "NOTCH3"
samples <- c(liver, kidney)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H7/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}
#PDGFRB
gene_name <- "PDGFRB"
samples <- c(bladder)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H7/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}
### H9
#STAT3
gene_name <- "STAT3"
samples <- c(prostate, bladder)
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  sce$gene <- assay(sce)[gene_name,]
  
  v <- .make_triangle_subspots(colData(sce), fill = "gene")
  p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
    scale_fill_viridis_c(option = "B") + labs(fill = "")
  pdf(paste0("Desktop/markers/H9/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
  plot(p)
  dev.off()
  rm(sce)
  gc()
}


seurat <- readRDS("Desktop/enhanced/new/pdac/P306_enhanced.rds")

ST_sce <- SingleCellExperiment(assays =list(SCT = as.matrix(GetAssayData(seurat, slot = "data"))), colData=seurat@meta.data)
metadata(ST_sce)$BayesSpace.data <- list()
metadata(ST_sce)$BayesSpace.data$platform <- "Visium"
metadata(ST_sce)$BayesSpace.data$is.enhanced <- TRUE
saveRDS(ST_sce, "Desktop/enhanced/new/all/P306_enhanced.rds")



library(Data)
DFrame()

genes <- c("CDK1", "EGFR", "EZH2", "GAPDH", "POT1", "SUZ12", "TXN", "VCAN", "XRCC4")
for (sample in pdac) {
  print(sample)
  sce <- readRDS(paste0("Desktop/enhanced/new/all/", sample, "_enhanced.rds"))
  sce <- sce[,-1]
  for (gene_name in genes) {
    sce$gene <- assay(sce)[gene_name,]
    v <- .make_triangle_subspots(colData(sce), fill = "gene")
    p <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
      scale_fill_viridis_c(option = "B") + labs(fill = "")
    pdf(paste0("Desktop/markers/pdac/", gene_name, "/",  sample, "_", gene_name, ".pdf"))
    plot(p)
    dev.off()   
  }
  rm(sce)
  gc()
}


library(ggpubr)
ggboxplot(df_all[df_all$sample=="Colorectal",], y = "H11", x = "clusters", fill = "clusters")


sub <- df_all[df_all$sample=="Colorectal",]
library(data.table)
long <- melt(setDT(sub), measure.vars = paste0("H", 1:13), variable.name = "hallmark")
long$clusters <- as.factor(long$clusters)
ggboxplot(long, y = "value", x = "hallmark",fill = "clusters", palette = c("#1F78B4", "#B2DF8A", "#FF7F00", "#B15928", "#FFFF99"),
          xlab = "Hallmark", ylab = "Scaled hallmark activity") + scale_x_discrete(labels=hallmark_names) + theme(axis.text.x = element_text(angle = 315, hjust = 0),
                                                                                                                  plot.margin = margin(1,5,1,1, "cm"))
 