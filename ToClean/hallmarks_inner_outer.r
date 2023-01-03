STobject <- readRDS("Desktop/enhanced/Breast_enhanced.rds")

sd_val <- 1


factor1_spots <- STobject$factor_1 > mean(STobject$factor_1) +sd_val*sd(STobject$factor_1)

factor2_spots <- STobject$factor_2 > mean(STobject$factor_2) +sd_val*sd(STobject$factor_2)

factor3_spots <- STobject$factor_3 > mean(STobject$factor_3) +sd_val*sd(STobject$factor_3)

factor4_spots <- STobject$factor_4 > mean(STobject$factor_4) +sd_val*sd(STobject$factor_4)

factor5_spots <- STobject$factor_5 > mean(STobject$factor_5) +sd_val*sd(STobject$factor_5)

Tumor <- colnames(STobject)[factor5_spots]

Tumor <- c(colnames(STobject)[factor1_spots], colnames(STobject)[factor4_spots])
Tumor <- Tumor[!duplicated(Tumor)]

TME <- colnames(STobject)[factor5_spots]
TME <- c(colnames(STobject)[factor2_spots], colnames(STobject)[factor3_spots])
TME <- TME[!duplicated(TME)]

duplicated <- c(TME, Tumor)[duplicated(c(TME, Tumor))]

TME <- TME[!TME %in% duplicated]
Tumor <- Tumor[!Tumor %in% duplicated]

coord <- STobject@images[[1]]@coordinates
distances <- as.matrix(dist(coord, diag = T, upper = T))

distances <- distances[Tumor, TME]


hist(apply(distances,1,min))

min_dist <- apply(distances,1,min)

inner <- names(min_dist[min_dist > quantile(min_dist, 0.75)])
outer <- names(min_dist[min_dist < quantile(min_dist, 0.25)])


STobject$Cancer <- NA
STobject$Cancer[inner] <- "TME-distant"
STobject$Cancer[outer] <- "TME-proximal"

SpatialDimPlot(STobject, group.by = "Cancer", pt.size.factor = 0.65)

STobject@meta.data[,paste0("H",1:13)] <- scale(STobject@meta.data[,paste0("H",1:13)])

H_modules <- STobject@meta.data[,c("Cancer", paste0("H", 1:13))]
H_modules[2:14] <- scale(H_modules[2:14])
colnames(H_modules[2:14]) <- hallmark_names
mean(H_modules$H12)
sd(H_modules$H12)

library(reshape2)
H_modules_resh <- melt(H_modules,id.vars='Cancer', measure.vars=paste0("H", 1:13))
colnames(H_modules_resh)[1:2] <- c("Cancer regions", "Hallmark")
H_modules_resh <- H_modules_resh[!is.na(H_modules_resh$`Cancer regions`),]
H_modules_resh <- H_modules_resh[H_modules_resh$Hallmark %in% paste0("H", c(2,4,8,9,11,12)),]

H_modules_resh$Hallmark <- sapply(as.character(H_modules_resh$Hallmark), function(x) {
  n <- substr(x, 2, stop = nchar(as.character(x)))
  hallmark_names[as.numeric(n)]
})



library(ggpubr)
ggboxplot(H_modules_resh, x = "Cancer regions", y = "value", ylab = "Hallmark scores", merge = "flip",
                ggtheme = theme_classic(), fill = "Cancer regions", palette = "ucscgb",
                legend = "top",
                title = "",
                font.legend = c(15, "plain", "black"),
                font.tickslab = c(15,"black"),
                font.main = c(15, "plain", "black"),
                font.x = c(15, "plain", "black"),
                font.y = c(15, "plain", "black"), facet.by = "Hallmark"
) + stat_compare_means()

mean(H_modules_resh$value[H_modules_resh$cluster=="inner" & H_modules_resh$Hallmark == "H2"])-
mean(H_modules_resh$value[H_modules_resh$cluster=="outer" & H_modules_resh$Hallmark == "H2"])

#


min_dist <- apply(distances,2,min)

inner <- names(min_dist[min_dist >= quantile(min_dist, 0.75)])
outer <- names(min_dist[min_dist <= quantile(min_dist, 0.25)])


STobject$TME <- NA
STobject$TME[inner] <- "Cancer-distant"
STobject$TME[outer] <- "Cancer-proximal"

SpatialDimPlot(STobject, group.by = "TME", pt.size.factor = 0.65)

STobject@meta.data[,paste0("H",1:13)] <- scale(STobject@meta.data[,paste0("H",1:13)])

H_modules <- STobject@meta.data[,c("TME", paste0("H", 1:13))]
H_modules[2:14] <- scale(H_modules[2:14])
mean(H_modules$H12)
sd(H_modules$H12)

H_modules_resh <- melt(H_modules,id.vars='TME', measure.vars=paste0("H", 1:13))
colnames(H_modules_resh)[1:2] <- c("TME regions", "Hallmark")
H_modules_resh <- H_modules_resh[!is.na(H_modules_resh$`TME regions`),]
H_modules_resh <- H_modules_resh[H_modules_resh$Hallmark %in% paste0("H", c(1,3,5,6,7,10,13)),]
H_modules_resh$Hallmark <- sapply(as.character(H_modules_resh$Hallmark), function(x) {
  n <- substr(x, 2, stop = nchar(as.character(x)))
  hallmark_names[as.numeric(n)]
})

ggboxplot(H_modules_resh, x = "TME regions", y = "value", ylab = "Hallmark scores", merge = "flip",
          ggtheme = theme_classic(), fill = "TME regions", palette = "ucscgb",
          legend = "top",
          title = "",
          font.legend = c(15, "plain", "black"),
          font.tickslab = c(15,"black"),
          font.main = c(15, "plain", "black"),
          font.x = c(15, "plain", "black"),
          font.y = c(15, "plain", "black"), facet.by = "Hallmark"
) + stat_compare_means()

mean(H_modules_resh$value[H_modules_resh$cluster=="inner" & H_modules_resh$Hallmark == "H9"])
  mean(H_modules_resh$value[H_modules_resh$cluster=="outer" & H_modules_resh$Hallmark == "H9"])

  
  
  