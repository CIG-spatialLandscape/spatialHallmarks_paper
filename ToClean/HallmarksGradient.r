### Gradient method
library(Seurat)
library(ggplot2)
source("Desktop/IJC/datasets/IGTP/figuresPaper/scripts/utilities/plotfunct.r")
source("Desktop/IJC/datasets/IGTP/figuresPaper/scripts/utilities/AnnotateCorrelations.r")


df <- as.data.frame(matrix(nrow=0, ncol=5))

files <- list.files(path="Desktop/enhanced", pattern="*.rds", full.names=TRUE, recursive=FALSE)
files <- files[c(1,2,3,4,5,6,7,10,13,14,15,17,18,19)]

for (file in files) {
  STobject <- readRDS(file)
  print(file)
  
  
  name <- strsplit(strsplit(file, "/")[[1]][3], "_")[[1]][1] 
  if (length(strsplit(strsplit(file, "/")[[1]][3], "_")[[1]])==3) name <- paste0(name, "_FFPE")
  
  
  
  
  factor_spots <- lapply(1:5, function(x) {
    classification <- kmeans(scale(STobject@meta.data[,paste0("factor_", x)]), centers = 2)$cluster
    if (mean(STobject@meta.data[,paste0("factor_", x)][classification==1]) >
        mean(STobject@meta.data[,paste0("factor_", x)][classification==2])) return(colnames(STobject)[classification==1])
    else return(colnames(STobject)[classification==2])
  })
  
  
  cancer_factors <- which(annotation_factors[[name]] == "Cancer" | annotation_factors[[name]] == "Cancer-Stroma" )
  TME_factors <-which(annotation_factors[[name]] != "Cancer"  & annotation_factors[[name]] != "Cancer-Stroma" & annotation_factors[[name]] != "Unassigned" )

  cancer <- sapply(cancer_factors, function(x){
    factor_spots[[x]]
  })
  if (class(cancer) == "list") cancer <- do.call(c, cancer)
  cancer <- cancer[!duplicated(cancer)]
  
  TME <- sapply(TME_factors, function(x){
    factor_spots[[x]]
  })
  if (class(TME) == "list") TME <- do.call(c, TME)
  TME <- TME[!duplicated(TME)]
  
  #remove overlapped spots 
  duplicated <- c(TME, cancer)[duplicated(c(TME, cancer))]
  TME <- TME[!TME %in% duplicated]
  cancer <- cancer[!cancer %in% duplicated]
  
  
  #extract coordiantes
  coord <- STobject@images[[1]]@coordinates[c(cancer, TME),c("row", "col")]
  distances <- as.matrix(dist(coord, diag = T, upper = T))
  distances <- distances[cancer, TME]
  gc()  
  
  
  #Cancer
  min_dist <- apply(distances,1,min)
  
  for (hallmark in paste0("H", 1:13)) {
    tmp <- data.frame(d = min_dist, h=scale(STobject@meta.data[cancer,hallmark]))
    tmp$hallmark <- hallmark
    tmp$sample <- name
    tmp$h_type <- "Cancer"
    df <- rbind(df, tmp)
  }
  
  # 
  # for (hallmark in paste0("H", c(2,4,8,9,11,12))) {
  #   tmp <- data.frame(d = min_dist, h=scale(STobject@meta.data[cancer,hallmark]))
  #   tmp$hallmark <- hallmark
  #   tmp$sample <- name
  #   tmp$h_type <- "Cancer"
  #   df <- rbind(df, tmp)
  # }
  # 
  #TME
  min_dist <- apply(distances,2,min)
  for (hallmark in paste0("H", 1:13)) {
    tmp <- data.frame(d = min_dist, h=scale(STobject@meta.data[TME,hallmark]))
    tmp$hallmark <- hallmark
    tmp$sample <- name
    tmp$h_type <- "TME"
    df <- rbind(df, tmp)
  }
  # for (hallmark in paste0("H", c(1,3,5,6,7,10,13))) {
  #   tmp <- data.frame(d = min_dist, h=scale(STobject@meta.data[TME,hallmark]))
  #   tmp$hallmark <- hallmark
  #   tmp$sample <- name
  #   tmp$h_type <- "TME"
  #   df <- rbind(df, tmp)
  # }
  
  rm(STobject, coord, distances)
  gc()
}


head(df)

df <- read.table("Desktop/IJC/datasets/IGTP/figuresPaper/test.txt", sep = "\t", header = T)


sub <- df[df$sample=="Ductal_FFPE", ]
sub <- df[(df$sample=="Ductal_FFPE" | df$sample=="OV4A" | df$sample=="Ovarian") & df$hallmark=="H12", ]
sub <- df[df$hallmark %in% paste0("H", c(2,4,8,9,11,12)), ]
ggplot(sub, aes(x=d, y=h, col=sample)) + geom_smooth() +
  facet_grid(~h_type, scales = "free") + scale_color_manual(values = DiscretePalette(13))
ggplot(sub, aes(x=d, y=h, col=sample)) + geom_smooth(method = "loess") + facet_grid(~sample, scales = "free")

ggplot(sub, aes(x=d, y=h, col=sample)) + geom_point()
unique(sub$sample)

ggplot(sub, aes(x=d, y=h, col=sample)) + geom_smooth(method = "loess") + facet_grid(rows = vars(hallmark), cols=vars(sample), scales = "free")


STobject$dist <- min_dist
SpatialFeaturePlot(STobject, features = c("dist", "H12"))



sub <- df[df$hallmark %in% paste0("H", c(1,3,5,6,7,10,13)), ]
sub <- df[df$h_type == "TME", ]

ggplot(sub, aes(x=d, y=h, col=sample)) + geom_smooth(method = "loess") + facet_grid(rows = vars(hallmark), cols=vars(sample), scales = "free")
