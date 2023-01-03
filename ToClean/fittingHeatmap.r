STobject <- readRDS("Desktop/enhanced/toUpload/Ovarian_enhanced.rds")

source("Desktop/IJC/datasets/IGTP/figuresPaper/scripts/utilities/plotfunct.r")

mat <- read.table("Desktop/distance/raw/matrix.txt", sep = "\t")


mat <- mat[mat$sample == "Ovarian" & mat$hallmark == "H1",]

mat$scaled_d <- NA
for (sample in unique(mat$sample)) {
  for (compartment in unique(mat$h_type)) {
    mat$scaled_d[mat$sample == sample & mat$h_type == compartment] <- scales::rescale(mat$d[mat$sample == sample & mat$h_type == compartment], to=c(0,1))
  }
}

mat$scaled_d[mat$h_type=="TME"] <- -mat$scaled_d[mat$h_type=="TME"]

rows <- sapply(rownames(mat), function(x) {
  s <- strsplit(x, split = ".", fixed = T)[[1]]
  paste0(s[1], ".", substr(s[2], start = 1, stop = 1))
})


names(rows) <- NULL

STobject$distances <- NA
STobject@meta.data[rows, "distances"] <- mat$scaled_d

SpatialFeaturePlot(STobject, features = c("distances", "H2"), pt.size.factor = 0.65)

STobject@meta.data[,paste0("H", 1:13)] <- scale(STobject@meta.data[,paste0("H", 1:13)])

STobject$score <- abs(STobject$distances)*STobject$H5
SpatialFeaturePlot(STobject, features = c("distances", "H5","score"), pt.size.factor = 0.65)





#############

mat <- read.table("Desktop/distance/all/new_distances.txt", sep = "\t")

test <- mat[mat$hallmark == "H2",]
test$scaled_d[test$h_type=="TME"] <- -test$scaled_d[test$h_type=="TME"]
test$d[test$h_type=="TME"] <- -test$d[test$h_type=="TME"]


test <- mat[mat$hallmark == "H12" & mat$h_type %in% c("Cancer_TME", "TME_Cancer") & mat$d <= 500,]
test$scaled_d[test$h_type=="TME_Cancer"] <- -test$scaled_d[test$h_type=="TME_Cancer"]
test$d[test$h_type=="TME_Cancer"] <- -test$d[test$h_type=="TME_Cancer"]


fitting2 <- as.data.frame(matrix(nrow=0, ncol=7))
for (tumor in unique(test$sample)) {
  print(tumor)
  for (hallmark in unique(test$hallmark)) {
    test1 <- test[test$sample==tumor & test$hallmark==hallmark,]
    model <- loess(h ~ d, data=test1)
    xseq <- seq(from=-500, 500, by=1)
    pred <- predict(model, newdata = data.frame(d = xseq), se=T)
    y = pred$fit
    ci <- pred$se.fit * qt(0.95 / 2 + .5, pred$df)
    ymin = y - ci
    ymax = y + ci
    loess.test <- data.frame(x = xseq, y, ymin, ymax, se = pred$se.fit)
    loess.test$sample <- tumor
    loess.test$hallmark <- hallmark
    fitting2 <- rbind(fitting2, loess.test)
  }
}
