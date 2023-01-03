source()
plotHIRES <- function(sample, feature,
                      compartment="whole", radar=FALSE,
                      palette=NULL, both=FALSE,
                      reverse_axis = FALSE,
                      shift_x = 0,
                      shift_y = 0
                      ) {
  sce <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/", sample, "_sce_enhanced.rds"))
  sce <- sce[,-1]
  if (!radar) {
    df <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/output_df/", sample, ".txt"))
    sce[, feature] <- df[ ,feature]
    if (compartment=="TME") clusters <- 4:5
    else if (compartment=="Cancer") clusters <- 1:2
    else clusters <- 1:5
    v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% clusters,], fill = feature)
    ref_v <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/", sample, ".rds"))
    hires <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/", sample, ".rds"))
    for (spot in unique(v$spot)) {
      v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
      v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
    }
    if (is.null(palette)) palette <- rev(paletteer_c("grDevices::YlOrRd", 30))
    
    plot <- hires + geom_polygon(data=v,  aes_(x=~imagerow+shift_x, y=~imagecol+shift_y, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
      scale_fill_gradientn(feature, colours = palette) 
    return(plot)
  }
  
}


sample <- "CUP295"
sce <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/CUP295_sce_enhanced.rds")
sce <- sce[,-1]
sce$H2 <- df_all[df_all$sample==sample, "H11"]
sce$H7 <- df_all[df_all$sample==sample, "H6"]
sce$estimate <- df_all[df_all$sample==sample, "estimate"]
sce$clusters <- df_all[df_all$sample==sample, "clusters"]


v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill = "H2")
v2 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(4,5),], fill = "H7")

ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/Ovarian_FFPE.rds")
hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/Ovarian_FFPE.rds")
for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}


for (spot in unique(v2$spot)) {
  v2$imagecol[v2$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v2$imagerow[v2$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}



hires + geom_polygon(data=v2,  aes_(x=~imagecol, y=~imagerow, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H3", colours = rev(paletteer_c("grDevices::YlOrRd", 30))  ) + 
  ggnewscale::new_scale("fill") + geom_polygon(data=v,  aes_(x=~imagecol, y=~imagerow + 60, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H2", colours = rev(
    
    paletteer_c("grDevices::Blues 2", 30)
  ))