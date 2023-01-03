library(ggpointdensityplot)


files <- list.files("Desktop/IJC/datasets/IGTP/figuresPaper/distances_fixed/", full.names = T)
mat <- lapply(files, function(file){
  read.table(file, sep = "\t", header = T)
})
mat <- do.call(rbind, mat)

mat <- mat[mat$d <= 500, ]
mat$scaled_d <- NA
mat$scaled_h <- NA
for (sample in unique(mat$sample)) {
  for (compartment in unique(mat$h_type)) {
    #mat$scaled_d[mat$sample == sample & mat$h_type == compartment] <- scales::rescale(mat$d[mat$sample == sample & mat$h_type == compartment], to=c(0.001,1))
    mat$scaled_d[mat$sample == sample & mat$h_type == compartment] <- scales::rescale(mat$d[mat$sample == sample & mat$h_type == compartment], to=c(0.001,1))
    mat$scaled_h[mat$sample == sample & mat$h_type == compartment] <- scale(mat$h[mat$sample == sample & mat$h_type == compartment])
    }
}
for (sample in unique(mat$sample)) {
  for (compartment in unique(mat$h_type)) {
    #mat$scaled_d[mat$sample == sample & mat$h_type == compartment] <- scales::rescale(mat$d[mat$sample == sample & mat$h_type == compartment], to=c(0.001,1))
    mat$scaled_d[mat$sample == sample & mat$h_type == compartment] <- scales::rescale(mat$d[mat$sample == sample & mat$h_type == compartment], to=c(0.001,1))
    mat$scaled_h[mat$sample == sample & mat$h_type == compartment] <- scale(mat$h[mat$sample == sample & mat$h_type == compartment])
  }
}


test <- mat[mat$sample == "OV4A" & mat$hallmark == "H12" & mat$h_type=="Cancer_TME",]
test$h <- scale(test$h)
ggplot(test, aes(x=d, y = h)) + geom_point() + geom_smooth(method = "lm") 

test <- mat[mat$sample == "Colorectal" & mat$hallmark == "H12" & mat$h_type=="Cancer_TME",]
test$h <- scale(test$h)
ggplot(test, aes(x=d, y = h)) + geom_point() + geom_smooth(method = "lm") 

test <- mat[mat$sample == "DU2" & mat$hallmark == "H12" & mat$h_type=="Cancer_TME",]
test$h <- scale(test$h)
ggplot(test, aes(x=d, y = h)) + geom_point() + geom_smooth(method = "lm") 


test <- mat[mat$hallmark == "H12",]
test$d[test$h_type=="TME"] <- -test$d[test$h_type=="TME"]
ggplot(test, aes(x=d, y = h)) + geom_point() + geom_smooth() + facet_wrap(~hallmark)


test$scaled_d[test$h_type=="TME"] <- -test$scaled_d[test$h_type=="TME"]
test$d[test$h_type=="TME"] <- -test$d[test$h_type=="TME"]
ggplot(test[test$scaled_d < 0.75 & test$scaled_d > -0.75,], aes(x=scaled_d, y = h)) + 
  geom_smooth() + geom_vline(xintercept = 0, col="black", linetype='dashed')+ geom_smooth() + theme_classic() + facet_wrap(~hallmark)

ggplot(test[test$scaled_d < 0.75 & test$scaled_d > -0.75,], aes(x=d, y = h)) + geom_smooth() 
  
  ggplot(test[test$scaled_d < 0.75 & test$scaled_d > -0.75,], aes(x=scaled_d, y = h)) + geom_vline(xintercept = 0, col="black", linetype='dashed')+ geom_smooth() + theme_classic()

  library(mgcv)
  

  summary(mod_gam1) 
  
  
  
  
  
mat <- df_hpa
source("Desktop/IJC/datasets/IGTP/figuresPaper/scripts/utilities/AnnotateCorrelations.r")
for (hallmark in 1:13) {
  mat[mat$hallmark == paste0("H", hallmark), "hallmark"] <- hallmark_names[hallmark]
}
#mat <- read.table("Desktop/IJC/datasets/IGTP/figuresPaper/distances_fixed/Colorectal_distances.txt", sep = "\t", header = T)
hallmark <- 12
threshold <- 10000
for (hallmark in 1:13) {
  test1 <- mat[mat$h_type %in% c("TME_Cancer", "Cancer_TME")  & mat$d <= threshold & mat$hallmark == hallmark_names[hallmark],]
  test1$d[test1$h_type=="TME_Cancer"] <- -test1$d[test1$h_type=="TME_Cancer"]
  test1$group <- "TME to Cancer"
  mod_gam1 = gam(h ~ s(d, bs = "cs"), data = test1, method = "REML", )
  df1 <- data.frame(d=seq(min(test1$d), max(test1$d), by=0.1))
  df1$prediction <- predict(mod_gam1, newdata = df1) 
  df1$compartment <- "TME"
  df1$compartment[df1$d > 0] <- "Cancer"
  df1$group <- "TME to Cancer"
  
  test2 <- mat[mat$h_type %in% c("TME_Buffer", "Buffer_TME")  & mat$d <= threshold & mat$hallmark == hallmark_names[hallmark],]
  test2$d[test2$h_type=="Buffer_TME"] <- -test2$d[test2$h_type=="Buffer_TME"]
  test2$group <- "Buffer to TME"
  mod_gam2 = gam(h ~ s(d, bs = "cs"), data = test2, method = "REML", )
  df2 <- data.frame(d=seq(min(test2$d), max(test2$d), by=0.1))
  df2$prediction <- predict(mod_gam2, newdata = df2) 
  df2$compartment <- "Buffer"
  df2$compartment[df2$d > 0] <- "TME"
  df2$group <- "Buffer to TME"
  
  test3 <- mat[mat$h_type %in% c("Cancer_Buffer", "Buffer_Cancer")  & mat$d <= threshold & mat$hallmark == hallmark_names[hallmark],]
  test3$d[test3$h_type=="Cancer_Buffer"] <- -test3$d[test3$h_type=="Cancer_Buffer"]
  test3$group <- "Cancer to Buffer"
  mod_gam3 = gam(h ~ s(d, bs = "cs"), data = test3, method = "REML", )
  df3 <- data.frame(d=seq(min(test3$d), max(test3$d), by=0.1))
  df3$prediction <- predict(mod_gam3, newdata = df3) 
  df3$compartment <- "Cancer"
  df3$compartment[df3$d > 0] <- "Buffer"
  df3$group <- "Cancer to Buffer"
  
  #test <- rbind(test1, test2, test3)
  #test$group <- factor(test$group, levels = c("Cancer to Buffer", "Buffer to TME", "TME to Cancer"))
  df <- rbind(df1, df2, df3)
  df$group <- factor(df$group, levels = c("Cancer to Buffer", "Buffer to TME", "TME to Cancer"))
  p <- ggplot(df, aes(x=d, y=prediction, col = compartment))  + geom_line(size=2) + 
    geom_vline(xintercept = 0, col="black", linetype='dashed')  + #+ scale_x_continuous(n.breaks = 3, labels = c(1000,  0,  1000))
    theme_classic() + labs(x ="Distance (microns)", y = "Hallmark activity")  + facet_wrap(~group, strip.position = "bottom") + ggtitle(hallmark_names[hallmark]) + 
    theme(plot.title = element_text(hjust = 0.5, size = 25), axis.text = element_text(size=15), axis.text.x = element_text(angle = 315, hjust = 0),
          axis.title = element_text(size=30), legend.position = "none",
          strip.text.x = element_text(size = 15), plot.margin = margin(0.2,1,0,0, "cm")) +  geom_hline(yintercept = 0, linetype='dashed') + 
    scale_color_manual(values = c(Cancer="orchid3",`Buffer` ="lightpink2", `TME`="lightgoldenrod2"))
  pdf(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/fittings/H", hallmark, ".pdf"), height = 6)
  plot(p)
  dev.off()  
}  
    

  test1 <- mat[mat$h_type %in% c("TME_Cancer", "Cancer_TME")  & mat$d <= 500 & mat$hallmark == hallmark_names[hallmark],]
  test1$d[test1$h_type=="TME_Cancer"] <- -test1$d[test1$h_type=="TME_Cancer"]
mod_gam1 = gam(h ~ s(d, bs = "cs"), data = test1, method = "REML", )
df <- data.frame(d=seq(-500, 500, by=0.1))
df$prediction <- predict(mod_gam1, newdata = df) 
df$group <- "1"
df$group[df$d > 0] <- "2"
ggplot(df, aes(x=d, y=prediction, col=group)) + geom_line()  + 
  geom_vline(xintercept = 0, col="black", linetype='dashed')  +
  theme_classic() + labs(x ="Distance (μm)", y = "Hallmark activity") 




test <- mat[mat$h_type %in% c("Cancer_Buffer", "Buffer_Cancer")  & mat$d <= 500,]
test$d[test$h_type=="Buffer_Cancer"] <- -test$d[test$h_type=="Buffer_Cancer"]
ggplot(test, aes(x=d, y=h))  + geom_smooth() + facet_wrap(~hallmark, scales="free")


  
  
  ################################################################################
#                             Distance map
################################################################################


STobject <- readRDS("Desktop/enhanced/toUpload/Colorectal_enhanced.rds")

mat <- read.table("Desktop/distance/all/new_distances.txt", sep = "\t")


test <- mat[mat$sample == "Colorectal" &  & mat$hallmark == "H12",]

#mat$scaled_d <- NA
#for (sample in unique(mat$sample)) {
#  for (compartment in unique(mat$h_type)) {
#    mat$scaled_d[mat$sample == sample & mat$h_type == compartment] <- scales::rescale(mat$d[mat$sample == sample & mat$h_type == compartment], to=c(0,1))
#    }
#}

test$scaled_d[test$h_type=="TME"] <- -test$scaled_d[test$h_type=="TME"]

rows <- sapply(rownames(test), function(x) {
  s <- strsplit(x, split = ".", fixed = T)[[1]]
  paste0(s[1], ".", substr(s[2], start = 1, stop = 1))
})


names(rows) <- NULL

STobject$distances <- NA
STobject@meta.data[rows, "distances"] <- test$scaled_d
library(Seurat)
SpatialFeaturePlot(STobject, features = "distances", pt.size.factor = 0.65)
SpatialDimPlot(STobject, pt.size.factor = 0.65)

library(paletteer)
source("Desktop/IJC/datasets/IGTP/figuresPaper/scripts/utilities/plotfunct.r")
SpatialFeaturePlot(STobject, features = "distances", pt.size.factor = 0.6) + scale_fill_gradientn("", colours = rev(paletteer_c("ggthemes::Orange-Blue-White Diverging", 30) ),breaks=c(-1,0,1), labels=c("TME core", "Interface","Cancer core")) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),legend.key.height = unit(3, 'cm'), 
        legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25), legend.margin=margin(l=-10))  + SpatialFeaturePlot(STobject, features = "H12", pt.size.factor = 0.6) + scale_fill_gradientn("", colours = viridis::inferno(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),legend.key.height = unit(3, 'cm'), 
        legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25), legend.margin=margin(l=-10))  + ggtitle(hallmark_names[12])

source()
STobject$H12 <- scale(STobject$H12)
SpatialFeaturePlot(STobject, features = "H12", pt.size.factor = 0.6) + scale_fill_gradientn("", colours = viridis::inferno(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),legend.key.height = unit(3, 'cm'), 
        legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25), legend.margin=margin(l=-10))  + ggtitle(hallmark_names[12])


