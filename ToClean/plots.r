
source("Desktop/IJC/datasets/IGTP/figuresPaper/scripts/utilities/BayesSpace_functions.r")
library(paletteer)


samples <- c("P3","P4","P6","P7","P8",paste0("M", 1:4))
palette <- RColorBrewer::brewer.pal(12, "Paired")
for (sample in samples) {
  sce <- readRDS(paste0("Desktop/IJC/datasets/Public/", sample, "/RDS/", sample, "_sce.rds"))
  clusterPlot(sce, palette = palette) + coord_fixed() + coord_flip() + scale_x_reverse() + scale_y_reverse() + theme(legend.position = "none")
  ggsave(paste0("Desktop/plots/", sample, "_cluster.png"), bg = "white", width = 7, height = 7)
  rm(sce)
  gc()
}


head(df_all)
sub <- df_all[df_all$sample=="OV4A",]
plot(sub$score, sub$estimate)
cor(sub$score, sub$estimate)

c <- c()
for (sample in unique(df_all$sample)) {
  sub <- df_all[df_all$sample==sample,]
  c <- rbind(c, c(sample, cor(sub$score, sub$estimate)))
}



sample <- "Co2"

#df_tme <- read.table("Desktop/df_tme58.txt", sep = "\t")
sce <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/Co2_sce_enhanced.rds")
sce <- sce[,-1]
sce$cnv <- df_all[df_all$sample==sample, "cnv_cluster"]
v <- .make_triangle_subspots(colData(sce), fill = "cnv")
ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~as.factor(fill))) +  theme_void() + coord_equal() 
#+ coord_flip() + scale_y_reverse()



v <- .make_triangle_subspots(colData(sce), fill = "estimate")
ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  scale_fill_gradientn("", colours = rev(paletteer_c("ggthemes::Orange-Blue Diverging", 30) ))

v <- .make_triangle_subspots(colData(sce), fill = "score")
ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("", colours = rev(paletteer_c("ggthemes::Orange-Blue Diverging", 30) ))

v <- .make_triangle_subspots(colData(sce), fill = "H11")
ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("", colours = rev(paletteer_c("ggthemes::Orange-Blue Diverging", 30) ))

v <- .make_triangle_subspots(colData(sce), fill = "clusters")
ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~factor(fill))) +  theme_void() + coord_equal()
clusterPlot(sce)


######3

v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill = "H11")

sub  <- colData(sce)[sce$clusters %in% c(4,5),]
sub$H11_cancer <- df_tme[df_tme$sample==sample, "H11_cancer"]
sub$H8_cancer <- df_tme[df_tme$sample==sample, "H8_cancer"]
sub$H2_cancer <- df_tme[df_tme$sample==sample, "H2_cancer"]

v2 <- .make_triangle_subspots(sub, fill = "H11_cancer")



p1 <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("", colours = rev(paletteer_c("ggthemes::Orange-Blue Diverging", 30) )) 


v <- .make_triangle_subspots(colData(sce), fill = "H11")
ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("", colours = rev(paletteer_c("ggthemes::Orange-Blue Diverging", 30) ))

v2 <- .make_triangle_subspots(sub, fill = "H11_cancer")
p1 <- ggplot()  + geom_polygon(data=v2,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("", colours = rev(paletteer_c("ggthemes::Orange-Blue Diverging", 30) ))


v2 <- .make_triangle_subspots(sub, fill = "H2_cancer")
p1 <- ggplot()  + geom_polygon(data=v2,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("", colours = rev(paletteer_c("ggthemes::Orange-Blue Diverging", 30) ))


v2 <- .make_triangle_subspots(sub, fill = "H8_cancer")
 ggplot()  + geom_polygon(data=v2,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("", colours = rev(paletteer_c("ggthemes::Orange-Blue Diverging", 30) ))



v2 <- .make_triangle_subspots(sub, fill = "H3")
p2 <- ggplot()  + geom_polygon(data=v2,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("", colours = rev(paletteer_c("ggthemes::Orange-Blue Diverging", 30) ))
p1+p2


###
v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill = "H2")
v2 <- .make_triangle_subspots(sub, fill = "H3")

new_scale <- function(new_aes) {
  structure(ggplot2::standardise_aes_names(new_aes), class = "new_aes")
}


for (spot in unique(v$spot)) {
    v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
    v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}


for (spot in unique(v2$spot)) {
  v2$imagecol[v2$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v2$imagerow[v2$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}

hires  + geom_polygon(data=v2,  aes_(x=~imagecol, y=~imagerow, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H3", colours = rev(paletteer_c("grDevices::YlOrRd", 30))  ) + 
  new_scale("fill") + geom_polygon(data=v,  aes_(x=~imagecol, y=~imagerow, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H2", colours = rev(
    
    paletteer_c("grDevices::Blues 2", 30)
  ))




##########
sample <- "M3"
sce <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/M3_sce_enhanced.rds")
sce <- sce[,-1]
sce$H2 <- df_all[df_all$sample==sample, "H2"]
sce$H12 <- df_all[df_all$sample==sample, "H12"]
sce$H11 <- df_all[df_all$sample==sample, "H11"]
sce$H9 <- df_all[df_all$sample==sample, "H9"]
sce$estimate <- df_all[df_all$sample==sample, "estimate"]
sce$clusters <- df_all[df_all$sample==sample, "clusters"]

ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/M3.rds")
hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/M3.rds")
v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill = "H11")
v2 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(4,5),], fill = "H9")

for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}


for (spot in unique(v2$spot)) {
  v2$imagecol[v2$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v2$imagerow[v2$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}


hires + geom_polygon(data=v2,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H9", colours = rev(paletteer_c("grDevices::YlOrRd", 30))  ) + 
  ggnewscale::new_scale("fill") + geom_polygon(data=v,  aes_(x=~imagerow+20, y=~imagecol + 100, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H12", colours = rev(
    
    paletteer_c("grDevices::Blues 2", 30)
  ))

ggsave("Desktop/P3_H12H9.png", bg = "white", width = 7, height = 7)

d1 <- v[v$spot %in% rownames(colData(sce)[sce$clusters %in% c(1,2),]),]
for (spot in unique(d1$spot)) {
  d1$fill[d1$spot == spot] <- colData(sce)[spot, "H11"] 

}
d2 <- v[v$spot %in% rownames(colData(sce)[sce$clusters %in% c(4,5),]),]
ggplot()  + geom_polygon(data=d2,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H13", colours = rev(paletteer_c("grDevices::YlOrRd", 30))  ) + 
  new_scale("fill") + geom_polygon(data=d1,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H11", colours = rev(
    
    paletteer_c("grDevices::Blues 2", 30)
  ))


##########
sample <- "P3"
sce <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/P3_sce_enhanced.rds")
sce <- sce[,-1]
sce$H2 <- df_all[df_all$sample==sample, "H2"]
sce$H11 <- df_all[df_all$sample==sample, "H11"]
sce$H13 <- df_all[df_all$sample==sample, "H13"]
sce$estimate <- df_all[df_all$sample==sample, "estimate"]
sce$clusters <- df_all[df_all$sample==sample, "clusters"]


v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill = "H11")
v2 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(4,5),], fill = "H13")

for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}


for (spot in unique(v2$spot)) {
  v2$imagecol[v2$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v2$imagerow[v2$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}


hires + geom_polygon(data=v2,  aes_(x=~imagerow+20, y=~imagecol+60, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H13", colours = rev(paletteer_c("grDevices::YlOrRd", 30))  ) + 
  ggnewscale::new_scale("fill") + geom_polygon(data=v,  aes_(x=~imagerow+20, y=~imagecol + 60, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H11", colours = rev(
    
   paletteer_c("grDevices::Blues 2", 30)
  ))

ggsave("Desktop/P3_H11H13.png", bg = "white", width = 7, height = 7)


ggplot()  + geom_polygon(data=v2,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H13", colours = rev(paletteer_c("grDevices::YlOrRd", 30))  ) + 
  new_scale("fill") + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H11", colours = rev(
    
    paletteer_c("grDevices::Blues 2", 30)
  ))

##
##########
sample <- "Co3"
sce <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/Co3_sce_enhanced.rds")
sce <- sce[,-1]
sce$H2 <- df_all[df_all$sample==sample, "H2"]
sce$H7 <- df_all[df_all$sample==sample, "H7"]
sce$H13 <- df_all[df_all$sample==sample, "H13"]
sce$estimate <- df_all[df_all$sample==sample, "estimate"]
sce$clusters <- df_all[df_all$sample==sample, "clusters"]


v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill = "H2")
v2 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(4,5),], fill = "H7")

for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}


for (spot in unique(v2$spot)) {
  v2$imagecol[v2$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v2$imagerow[v2$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}
ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/Co3.rds")
hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/Co3.rds")

hires + geom_polygon(data=v2,  aes_(x=~imagerow, y=~imagecol+100, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H3", colours = rev(paletteer_c("grDevices::YlOrRd", 30))  ) + 
  ggnewscale::new_scale("fill") + geom_polygon(data=v,  aes_(x=~imagerow+20, y=~imagecol + 60, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H2", colours = rev(
    
    paletteer_c("grDevices::Blues 2", 30)
  ))

ggsave("Desktop/P3_H11H13.png", bg = "white", width = 7, height = 7)


ggplot()  + geom_polygon(data=v2,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H13", colours = rev(paletteer_c("grDevices::YlOrRd", 30))  ) + 
  new_scale("fill") + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H11", colours = rev(
    
    paletteer_c("grDevices::Blues 2", 30)
  ))


##
##########
sample <- "OVFFPE"
sce <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/Ovarian_FFPE_sce_enhanced.rds")
sce <- sce[,-1]
sce$H2 <- df_all[df_all$sample==sample, "H2"]
sce$H7 <- df_all[df_all$sample==sample, "H6"]
sce$H13 <- df_all[df_all$sample==sample, "H13"]
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

  ggsave("Desktop/P3_H11H13.png", bg = "white", width = 7, height = 7)


ggplot()  + geom_polygon(data=v2,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H13", colours = rev(paletteer_c("grDevices::YlOrRd", 30))  ) + 
  new_scale("fill") + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H11", colours = rev(
    
    paletteer_c("grDevices::Blues 2", 30)
  ))

##
##########
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

ggsave("Desktop/P3_H11H13.png", bg = "white", width = 7, height = 7)


ggplot()  + geom_polygon(data=v2,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H6", colours = rev(paletteer_c("grDevices::YlOrRd", 30))  ) + 
  new_scale("fill") + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H11", colours = rev(
    
    paletteer_c("grDevices::Blues 2", 30)
  ))

##
##########
sample <- "DU13"
sce <- readRDS("Desktop/IJC/datasets/Dutreneo/RDS/enhancement/DU13_sce.rds")
sce <- sce[,-1]
sce$H2 <- df_all[df_all$sample==sample, "H11"]
sce$H7 <- df_all[df_all$sample==sample, "H1"]
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

ggsave("Desktop/P3_H11H13.png", bg = "white", width = 7, height = 7)


ggplot()  + geom_polygon(data=v2,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H13", colours = rev(paletteer_c("grDevices::YlOrRd", 30))  ) + 
  new_scale("fill") + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H11", colours = rev(
    
    paletteer_c("grDevices::Blues 2", 30)
  ))

##



##########
sample <- "cHC1T"
sce <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/cHC-1T_sce_enhanced.rds")
sce <- sce[,-1]
sce$H2 <- df_all[df_all$sample==sample, "H11"]
sce$H7 <- df_all[df_all$sample==sample, "H6"]
sce$H13 <- df_all[df_all$sample==sample, "H13"]
sce$estimate <- df_all[df_all$sample==sample, "estimate"]
sce$clusters <- df_all[df_all$sample==sample, "clusters"]


v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill = "H2")
v2 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(4,5),], fill = "H7")

ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/cHC-1T.rds")
hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/cHC-1T.rds")
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

ggsave("Desktop/P3_H11H13.png", bg = "white", width = 7, height = 7)


ggplot()  + geom_polygon(data=v2,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H13", colours = rev(paletteer_c("grDevices::YlOrRd", 30))  ) + 
  new_scale("fill") + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H11", colours = rev(
    
    paletteer_c("grDevices::Blues 2", 30)
  ))

###3


color_codes <-  list("Sustaining Proliferative Signaling"="#15CE59", 
                     "Evading Growth Suppressors" = "#701717",
                     "Avoiding Immune Destruction" = "#CB3BBD",
                     "Enabling Replicative Immortality" = "#4969D4",
                     "Tumour-Promoting Inflammation" = "#E6880D",
                     "Activating Invasion and Metastasis" = "#000000",
                     "Inducing Angiogenesis" = "#EE0F16",
                     "Genome Instability and Mutation" = "#132892",
                     "Resisting Cell Death" = "#8E909B",
                     "Deregulating Cellular Energetics" = "#71189E",
                     "Senescent cells" = "#05F3EB",
                     "Nonmutational Epigenetic reprogramming" = "#890269",
                     "Unlocking Phenotypic Plasticity" = "#95641A")
##########
sample <- "M3"
sce <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/M3_sce_enhanced.rds")
sce <- sce[,-1]
#Assign clusters based on NNMF activty
hallmarks <-  df_all[df_all$sample=="M3", paste0("H", 1:13)]
hallmarks <- scale(hallmarks)

#assign label to the maximum scaled factor activity
sce$h_clusters <- names(color_codes)[apply(hallmarks, 1, which.max)]
palette <- do.call(c, color_codes)


v <- .make_triangle_subspots(colData(sce), fill = "h_clusters")


ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/M3.rds")
hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/M3.rds")
for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}


v$fill <- factor(v$fill, levels = names(palette))
hires  + geom_polygon(data=v,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~fill)) +  
  theme_void() + scale_fill_manual(values = palette) + theme(legend.position = "none")

ggsave("Desktop/h_clusters_M3.png", bg = "white", width = 7, height = 7)


ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~as.factor(fill))) +  theme_void() + scale_fill_manual(values = palette)





##########
hallmark_names <-  list("H1" = "Sustaining Proliferative Signaling", 
                        "H2" = "Evading Growth Suppressors",
                        "H3" = "Avoiding Immune Destruction" ,
                        "H4" = "Enabling Replicative Immortality" ,
                        "H5" = "Tumour-Promoting Inflammation" ,
                        "H6" = "Activating Invasion and Metastasis",
                        "H7" = "Inducing Angiogenesis",
                        "H8" = "Genome Instability and Mutation",
                        "H9" = "Resisting Cell Death",
                        "H10" = "Deregulating Cellular Energetics",
                        "H11" = "Senescent cells" ,
                        "H12" = "Nonmutational Epigenetic reprogramming" ,
                        "H13" = "Unlocking Phenotypic Plasticity" )



sample <- "M3"
sce <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/M3_sce_enhanced.rds")
sce <- sce[,-1]
sce$clusters <- df_all[df_all$sample==sample, "clusters"]
sce <- sce[,sce$clusters %in% c(4,5)]

df_tme <- read.table("Desktop/df_tme58.txt", sep = "\t")

sce$H2_predictor <- df_tme[df_tme$sample==sample, "H2_cancer"]
sce$H4_predictor <- df_tme[df_tme$sample==sample, "H4_cancer"]
sce$H8_predictor <- df_tme[df_tme$sample==sample, "H8_cancer"]
sce$H11_predictor <- df_tme[df_tme$sample==sample, "H11_cancer"]
sce$H12_predictor <- df_tme[df_tme$sample==sample, "H12_cancer"]

sce$H1 <- df_all[df_all$sample==sample & df_all$clusters %in% c(4,5), "H1"]
sce$H3 <- df_all[df_all$sample==sample & df_all$clusters %in% c(4,5), "H3"]
sce$H5 <- df_all[df_all$sample==sample & df_all$clusters %in% c(4,5), "H5"]
sce$H6 <- df_all[df_all$sample==sample & df_all$clusters %in% c(4,5), "H6"]
sce$H7 <- df_all[df_all$sample==sample & df_all$clusters %in% c(4,5), "H7"]
sce$H9 <- df_all[df_all$sample==sample & df_all$clusters %in% c(4,5), "H9"]
sce$H13 <- df_all[df_all$sample==sample & df_all$clusters %in% c(4,5), "H13"]

ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/M3.rds")
hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/M3.rds")



for (h in paste0("H", c(1,3,5,6,7,9,13))) {
  print(h)
  v <- .make_triangle_subspots(colData(sce), fill = h)
  for (spot in unique(v$spot)) {
    v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
    v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
  }
  hires  + geom_polygon(data=v,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
    scale_fill_viridis_c(option = "B") + ggtitle(hallmark_names[[h]]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 25))
  ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RF_plots/", h, "_tme.png"), width = 6, height = 7)
}

for (h in paste0("H", c(2,4,8,11,12))) {
  print(h)
  v <- .make_triangle_subspots(colData(sce), fill = paste0(h, "_predictor"))
  for (spot in unique(v$spot)) {
    v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
    v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
  }
  hires  + geom_polygon(data=v,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
    scale_fill_viridis_c() + ggtitle(hallmark_names[[h]]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 25))
  ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RF_plots/", h, "_tme_predictor.png"), width = 6, height = 7)
}


sample <- "M3"
sce <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/M3_sce_enhanced.rds")
sce <- sce[,-1]
sce$clusters <- df_all[df_all$sample==sample, "clusters"]
sce <- sce[,sce$clusters %in% c(1,2)]

df_cancer2 <- read.table("Desktop/df_cancer.txt", sep = "\t")

sce$H1_predictor <- df_cancer[df_cancer$sample==sample, "H1_TME"]
sce$H3_predictor <- df_cancer[df_cancer$sample==sample, "H3_TME"]
sce$H5_predictor <- df_cancer[df_cancer$sample==sample, "H5_TME"]
sce$H6_predictor <- df_cancer[df_cancer$sample==sample, "H6_TME"]
sce$H7_predictor <- df_cancer[df_cancer$sample==sample, "H7_TME"]
sce$H9_predictor <- df_cancer[df_cancer$sample==sample, "H9_TME"]
sce$H13_predictor <- df_cancer[df_cancer$sample==sample, "H13_TME"]

sce$H2 <- df_all[df_all$sample==sample & df_all$clusters %in% c(1,2), "H2"]
sce$H4 <- df_all[df_all$sample==sample & df_all$clusters %in% c(1,2), "H4"]
sce$H8 <- df_all[df_all$sample==sample & df_all$clusters %in% c(1,2), "H8"]
sce$H11 <- df_all[df_all$sample==sample & df_all$clusters %in% c(1,2), "H11"]
sce$H12 <- df_all[df_all$sample==sample & df_all$clusters %in% c(1,2), "H12"]



ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/M3.rds")
hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/M3.rds")


for (h in paste0("H", c(2,4,8,11,12))) {
  print(h)
  v <- .make_triangle_subspots(colData(sce), fill = h)
  for (spot in unique(v$spot)) {
    v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
    v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
  }
  hires  + geom_polygon(data=v,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
    scale_fill_viridis_c(option = "B") + ggtitle(hallmark_names[[h]]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 25))
  ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RF_plots/", h, "_cancer_fix.png"), width = 6, height = 7)
}

for (h in paste0("H", c(1,3,5,6,7,9,13))) {
  print(h)
  v <- .make_triangle_subspots(colData(sce), fill = paste0(h, "_predictor"))
  for (spot in unique(v$spot)) {
    v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
    v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
  }
  hires  + geom_polygon(data=v,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
    scale_fill_viridis_c() + ggtitle(hallmark_names[[h]]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 25))
  ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RF_plots/", h, "_cancer_predictor_fix.png"), width = 6, height = 7)
}


hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/M3.rds")
v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/M3.rds")
v <- v[v$spot != "subspot_1.1",]
df_all$spotid <- sapply(rownames(df_all), function(spot){
  x <- strsplit(spot,  ".", fixed = T)
  paste0(x[[1]][1], ".", substr(x[[1]][2],1, 1))
})
sub_data <- df_all[df_all$sample == "M3",]
sub_data <- read.table("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/output_df/M3.txt")

library(paletteer)
d <- data.frame(estimate=sub_data$estimate, score=sub_data$score)
model <- lm(estimate~score, data = d)
sub_data$residuals <- model$residuals
library(paletteer)

cluster <- c("5" = "lightgoldenrod1", "4" = "lightgoldenrod3", "3" = "lightpink2", "2"= "orchid3", "1"= "orchid4")


v$fill <- sub_data[v$spot,"clusters"]
hires + geom_polygon(data=v,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~factor(fill))) +  theme_void() + coord_equal() +
  scale_fill_manual(values = rev(c("lightgoldenrod1", "lightgoldenrod3", "lightpink2",  "orchid3",  "orchid4"))) + labs(fill="") 
ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RF_plots/clusters.png"), width = 6, height = 7)

############

##########
sample <- "M3"
sce <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/M3_sce_enhanced.rds")
sce <- sce[,-1]
sce$clusters <- df_all[df_all$sample==sample, "clusters"]
sce$cnv <- cnv[cnv$sample == "M3", "cnv_cluster"]
v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill = "cnv")
ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~as.factor(fill))) +  theme_void() + coord_flip() + scale_x_reverse() + scale_y_reverse()

##
sample <- "M3"
sce <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/enhanced/",sample,"_sce_enhanced.rds"))
sce <- sce[,-1]
sce$clusters <- df_all[df_all$sample==sample, "clusters"]
sce$cnv <- cnv[cnv$sample == sample, "cnv_cluster"]
v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2) & sce$cnv != 3,], fill = "cnv")
ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~as.factor(fill))) +  theme_void() + coord_flip() + scale_x_reverse() + scale_y_reverse()
ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/M3.rds")
hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/M3.rds")
for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}
hires + geom_polygon(data=v,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~as.factor(fill))) +  theme_void() +
  labs(fill="CNV clones") +  scale_fill_manual(values = c("#E69F00", "#0072B2", "#009E73"))
####
sample <- "BreastA"
sce <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/enhanced/",sample,"_sce_enhanced.rds"))
sce <- sce[,-1]
sce$clusters <- df_all[df_all$sample==sample, "clusters"]
sce$cnv <- cnv[cnv$sample == sample, "cnv_cluster"]
v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2) & sce$cnv != 3,], fill = "cnv")
ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~as.factor(fill))) +  theme_void() + coord_flip() + scale_x_reverse() + scale_y_reverse()
ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/BreastA.rds")
hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/BreastA.rds")
for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}
hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~as.factor(fill))) +  theme_void() +
  labs(fill="CNV clones") +  scale_fill_manual(values = c("#E69F00", "#0072B2", "#009E73"))
####
sample <- "UKF260T"
sce <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/enhanced/",sample,"_sce_enhanced.rds"))
sce <- sce[,-1]
sce$clusters <- df_all[df_all$sample==sample, "clusters"]
sce$cnv <- df_all[df_all$sample == sample, "cnv_cluster"]
v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2) & sce$cnv != 3,], fill = "cnv")
ggplot()  + geom_polygon(data=v,  aes_(x=~-y.vertex, y=~-x.vertex, group=~spot, fill=~as.factor(fill))) +  theme_void() + coord_flip() + scale_x_reverse() + scale_y_reverse()
ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/UKF260T.rds")
hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/UKF260T.rds")
for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}
hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~as.factor(fill))) +  theme_void() +
  labs(fill="CNV clones") +  scale_fill_manual(values = c("#E69F00", "#0072B2", "#009E73"))

####
sample <- "HCC2T"
sce <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/enhanced/",sample,"_sce_enhanced.rds"))
sce <- sce[,-1]
sce$clusters <- df_all[df_all$sample==sample, "clusters"]
sce$cnv <- df_all[df_all$sample == sample, "cnv_cluster"]
v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2) & sce$cnv %in%  perc_filter$cnv_cluster[perc_filter$sample ==sample],], fill = "cnv")
ref_v <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/", sample, ".rds"))
hires <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/", sample, ".rds"))
for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}
pdf(paste0("Desktop/IJC/datasets/IGTP/figu"))
hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~as.factor(fill))) +  theme_void() +
  labs(fill="CNV clones") +  scale_fill_manual(values = c("#E69F00", "#0072B2", "#009E73"))
dev.off()
##########
sample <- "M3"
sce <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/M3_sce_enhanced.rds")
sce <- sce[,-1]
sce$H2 <- df_all[df_all$sample==sample, "H2"]
sce$H12 <- df_all[df_all$sample==sample, "H12"]
sce$H11 <- df_all[df_all$sample==sample, "H11"]
sce$H3 <- df_all[df_all$sample==sample, "H3"]
sce$estimate <- df_all[df_all$sample==sample, "estimate"]
sce$clusters <- df_all[df_all$sample==sample, "clusters"]

ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/M3.rds")
hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/M3.rds")
v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill = "H2")
v2 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(4,5),], fill = "H3")

for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}


for (spot in unique(v2$spot)) {
  v2$imagecol[v2$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v2$imagerow[v2$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}


hires + geom_polygon(data=v2,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H3", colours = rev(paletteer_c("grDevices::YlOrRd", 30))  ) + 
  ggnewscale::new_scale("fill") + geom_polygon(data=v,  aes_(x=~imagerow+20, y=~imagecol + 100, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn("H2", colours = rev(
    
    paletteer_c("grDevices::Blues 2", 30)
  ))

ggsave("Desktop/M3_H3H2.png", bg = "white", width = 7, height = 7)


###

sample <- "M2"
sce <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/M2_sce_enhanced.rds")
sce <- sce[,-1]
gtf <- read.table("Desktop/IJC/datasets/IGTP/figuresPaper/estimate/new_input/M2_exp.txt", sep = "\t")
markers <- c("ERBB2", "ESR1", "PGR", "MYC", "EPCAM", "KRT7", "KRT8", "KRT18")
markers <- c("ATM", "BARD1", "BRIP1", "CASP8", "CTLA4", "CYP19A1", "FGFR2", "H19", "LSP1", "MAP3K1", "MRE11A", "NBN", "RAD51", "TERT",
             "CD44", "ALDH1A1")

markers <- markers[markers %in% rownames(gtf)]
colData(sce) <- cbind(colData(sce), t(gtf[markers, rownames(colData(sce))]))

ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/M2.rds")
hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/M2.rds")

for (marker in markers) {
  v <- .make_triangle_subspots(colData(sce), fill = marker)
  
  for (spot in unique(v$spot)) {
    v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
    v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
  }
  
  hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
    scale_fill_viridis_c("",option = "B") + ggtitle(marker) + theme(plot.title = element_text(hjust = 0.5, size = 20))
  
  ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/M2_markers/M2_",marker, ".png"), bg = "white", width = 7, height = 7)  
}
colData(sce) <- cbind(colData(sce),  df_all[df_all$sample==sample, paste0("H",1:13)])
for (marker in paste0("H",1:13)) {
  v <- .make_triangle_subspots(colData(sce), fill = marker)
  
  for (spot in unique(v$spot)) {
    v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
    v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
  }
  
  hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
    scale_fill_viridis_c("",option = "B") + ggtitle(marker) + theme(plot.title = element_text(hjust = 0.5, size = 20))
  
  ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/M2_markers/M2_",marker, ".png"), bg = "white", width = 7, height = 7)  
}


###
sample <- "DU12"
sce <- readRDS("Desktop/IJC/datasets/Dutreneo/RDS/enhanced/DU12_sce_enhanced.rds")
sce <- sce[,-1]
gtf <- readRDS("Desktop/IJC/datasets/Dutreneo/RDS/enhancement/DU12_sce.rds")
markers <- c("ERBB2", "ESR1", "PGR", "MYC", "EPCAM", "KRT7", "KRT8", "KRT18")
colData(sce) <- cbind(colData(sce),  t(assays(gtf)[["SCT"]][markers,-1]))
colData(sce) <- cbind(colData(sce),  df_all[df_all$sample==sample, paste0("H",1:13)])
sce$estimate <- df_all[df_all$sample==sample, "estimate"] 
sce$clusters <- df_all[df_all$sample==sample, "clusters"]
sce$cnv <- cnv[cnv$sample == sample, "cnv_cluster"]
ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/DU12.rds")
hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/DU12.rds")

for (marker in markers) {
  v <- .make_triangle_subspots(colData(sce), fill = marker)
  
  for (spot in unique(v$spot)) {
    v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
    v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
  }
  
  hires + geom_polygon(data=v,  aes_(x=~imagerow+60, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
    scale_fill_viridis_c("",option = "B") + ggtitle(marker) + theme(plot.title = element_text(hjust = 0.5, size = 20))
  
  ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/DU12_markers/DU12_",marker, ".png"), bg = "white", width = 7, height = 7)  
}

for (marker in paste0("H",1:13)) {
  v <- .make_triangle_subspots(colData(sce), fill = marker)
  
  for (spot in unique(v$spot)) {
    v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
    v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
  }
  
  hires + geom_polygon(data=v,  aes_(x=~imagerow+60, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
    scale_fill_viridis_c("",option = "B") + ggtitle(marker) + theme(plot.title = element_text(hjust = 0.5, size = 20))
  
  ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/DU12_markers/DU12_",marker, ".png"), bg = "white", width = 7, height = 7)  
}

v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill = "cnv")
for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}
hires + geom_polygon(data=v,  aes_(x=~imagerow+60, y=~imagecol, group=~spot, fill=~as.factor(fill))) +  theme_void() +
  labs(fill="CNV clones")




###
sample <- "DU13"
sce <- readRDS("Desktop/IJC/datasets/Dutreneo/RDS/enhanced/DU13_sce_enhanced.rds")
sce <- sce[,-1]
gtf <- readRDS("Desktop/IJC/datasets/Dutreneo/RDS/enhancement/DU13_sce.rds")
markers <- c("ERBB2", "ESR1", "PGR", "MYC", "EPCAM", "KRT7", "KRT8", "KRT18")
colData(sce) <- cbind(colData(sce),  t(assays(gtf)[["SCT"]][markers,-1]))
colData(sce) <- cbind(colData(sce),  df_all[df_all$sample==sample, paste0("H",1:13)])
sce$estimate <- df_all[df_all$sample==sample, "estimate"] 
sce$clusters <- df_all[df_all$sample==sample, "clusters"]
sce$cnv <- cnv[cnv$sample == sample, "cnv_cluster"]
ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/DU12.rds")
hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/DU12.rds")

for (marker in markers) {
  v <- .make_triangle_subspots(colData(sce), fill = marker)
  
  
  ggplot() + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
    scale_fill_viridis_c("",option = "B") + ggtitle(marker) + theme(plot.title = element_text(hjust = 0.5, size = 20)) + coord_flip() + 
    scale_y_reverse() + scale_x_reverse()
  
  ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/DU13_markers/DU13_",marker, ".png"), bg = "white", width = 7, height = 7)  
}

for (marker in paste0("H",1:13)) {
  v <- .make_triangle_subspots(colData(sce), fill = marker)
  ggplot() + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
    scale_fill_viridis_c("",option = "B") + ggtitle(marker) + theme(plot.title = element_text(hjust = 0.5, size = 20))+  coord_flip() + 
    scale_y_reverse() + scale_x_reverse()
  
  ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/DU13_markers/DU13_",marker, ".png"), bg = "white", width = 7, height = 7)  
}


v <- .make_triangle_subspots(colData(sce), fill = "estimate")
ggplot() + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_viridis_c("",option = "B") + ggtitle(marker) + theme(plot.title = element_text(hjust = 0.5, size = 20))+  coord_flip() + 
  scale_y_reverse() + scale_x_reverse()

v <- .make_triangle_subspots(colData(sce), fill = "cnv")
ggplot() + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~as.factor(fill))) +  theme_void() +
  labs(fill="CNV clones")+  coord_flip() + 
  scale_y_reverse() + scale_x_reverse()

##

sample <- "M3"
sce <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/M3_sce_enhanced.rds")
sce <- sce[,-1]
sce$H10 <- df_all[df_all$sample==sample, "H10"]

ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/M3.rds")
hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/M3.rds")


  v <- .make_triangle_subspots(colData(sce), fill = "H10")
  
  for (spot in unique(v$spot)) {
    v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
    v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
  }
  
  hires + geom_polygon(data=v,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
    scale_fill_viridis_c("",option = "B") + ggtitle(marker) + theme(plot.title = element_text(hjust = 0.5, size = 20))
  
  ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/M3_markers/M3_",marker, ".png"), bg = "white", width = 7, height = 7)  

  
  ##
  
  sample <- "M1"
  sce <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/M1_sce_enhanced.rds")
  sce <- sce[,-1]
  sce$clusters <- as.factor(df_all[df_all$sample==sample, "clusters"])
  sce$comparments <- "Cancer"
  sce$comparments[sce$clusters %in% c(4,5)] <- "TME"
  sce$comparments[sce$clusters %in% c(3)] <- "Buffer"
  
  ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/M1.rds")
  hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/M1.rds")
  
  
  v <- .make_triangle_subspots(colData(sce), fill = "comparments")
  
  for (spot in unique(v$spot)) {
    v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
    v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
  }
  
  
  hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +
    theme_void() + coord_equal() +
    scale_fill_manual("Compartments", values = c("TME" = "lightgoldenrod1", "Buffer" = "lightpink2", "Cancer"= "orchid3"))
  
  ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/M3_markers/M3_",marker, ".png"), bg = "white", width = 7, height = 7)  
  
  
 ##
  
  sample <- "DU13"
  sce <- readRDS("Desktop/IJC/datasets/Dutreneo/RDS/enhancement/DU13_sce.rds")
  sce <- sce[,-1]
  sce$clusters <- as.factor(df_all[df_all$sample==sample, "clusters"])
  sce$comparments <- "Cancer"
  sce$comparments[sce$clusters %in% c(4,5)] <- "TME"
  sce$comparments[sce$clusters %in% c(3)] <- "Buffer"
  
  ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/DU13.rds")
  hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/DU13.rds")
  
  
  v <- .make_triangle_subspots(colData(sce), fill = "comparments")
  
  for (spot in unique(v$spot)) {
    v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
    v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
  }
  
  
  hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +
    theme_void() + coord_equal() +
    scale_fill_manual("Compartments", values = c("TME" = "lightgoldenrod1", "Buffer" = "lightpink2", "Cancer"= "orchid3"))
  
  ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/M3_markers/M3_",marker, ".png"), bg = "white", width = 7, height = 7)    ##
  
  
  #####
  
  sample <- "OVFFPE"
  sce <- readRDS(paste0("Desktop/enhanced/enhanced_sce/", sample, "_sce_enhanced.rds"))
  sce <- sce[,-1]
  sce$clusters <- df_all[df_all$sample==sample,  "clusters"]
  sce <- sce[,sce$clusters %in% c(4,5)]
  colData(sce)[,c(paste0("H",c(1,3,5,6,7,9,13)))] <- df_all[df_all$sample==sample & df_all$clusters %in% c(4,5), c(paste0("H",c(1,3,5,6,7,9,13)))]
  df_tme <- read.table("Desktop/df_tme_fix.txt")
  colData(sce)[,paste0("H", c(2,4,8,10,11,12), "_predictor")] <- df_tme[df_tme$sample==sample, paste0("H", c(2,4,8,10,11,12), "_cancer")]
  
  

  
  ref_v <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/",sample,".rds"))
  hires <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/",sample,".rds"))

  for (h in paste0("H", c(2,4,8,10,11,12))) {
    print(h)
    v <- .make_triangle_subspots(colData(sce), fill = paste0(h, "_predictor"))
    for (spot in unique(v$spot)) {
      v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
      v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
    }
    hires  + geom_polygon(data=v,  aes_(x=~imagecol, y=~imagerow, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
      scale_fill_viridis_c() + ggtitle(hallmark_names_list[[h]]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 25))
    ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RF_plots/OVFFPE/", h, "_tme_predictor.png"), width = 6, height = 7)
  }
  
  for (h in paste0("H", c(1,3,5,6,7,9,13))) {
    print(h)
    v <- .make_triangle_subspots(colData(sce), fill = h)
    for (spot in unique(v$spot)) {
      v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
      v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
    }
    hires  + geom_polygon(data=v,  aes_(x=~imagecol, y=~imagerow, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
      scale_fill_viridis_c(option = "B") + ggtitle(hallmark_names_list[[h]]) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 25))
    ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RF_plots/OVFFPE/", h, "_tme.png"), width = 6, height = 7)
  }
 
  
  
  
  sample <- "HCC1T"
  sce <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/HCC-1T_sce_enhanced.rds")
  sce <- sce[,-1]
  
  ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/HCC1T.rds")
  hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/HCC1T.rds")
  colData(sce) <- cbind(colData(sce),  df_all[df_all$sample==sample, c("clusters",paste0("H",1:13))])
  c(2,4,8,10,11,12)
  c(1,3,5,6,7,9,13)
  for (marker in paste0("H",c(1,3,5,6,7,9,13))) {
    v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% 4:5,], fill = marker)
    
    for (spot in unique(v$spot)) {
      v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
      v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
    }
    
    hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
      scale_fill_viridis_c("",option = "B") + ggtitle(marker) + theme(plot.title = element_text(hjust = 0.5, size = 20))
    
    ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/HCC1T_markers/HCC1T_",marker, ".png"), bg = "white", width = 7, height = 7)  
  }
