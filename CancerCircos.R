##################################################
## Project: Cancer Hallmarks
## Script purpose: Plot results of Random Forest models for Cancer Hallmarks
## Date: 22/12/2022
## Author: Sergi Cervilla & Mustafa Sibai
##################################################

library(circlize)
library(dplyr)
library(paletteer)
library(caret)
library(dplyr)
library(treeshap)
library(stringr)

############################## prediction of cancer compartment ###################

RF <- list()
#grab each Random Forest shapley model and put it in the same list
for(file in list.files("")) {
  print(file)
  RF[[strsplit(file, split = ".rds")[[1]]]] <- readRDS("")
}

### Extract Importances
RF_Imp.mt <- data.frame(matrix(nrow = 0, ncol = 18))

#extract importances of each RF model 
for (ID in names(RF)){
  #barplot of importances
  imp_plot <- plot_feature_importance(RF[[ID]], max_vars = 7)
  #function to scale from 0 to 100
  rescale <- function(x) (x-min(x))/(max(x) - min(x)) * 100
  #apply scaling to importance values
  imp_plot$data$importance <- rescale(imp_plot$data$importance)
  #set predictor as factor
  imp_plot$data$variable <- factor(imp_plot$data$variable, levels = str_sort(imp_plot$data$variable))
  #order by predictor
  imp_plot$data <- arrange(imp_plot$data, variable)
  #create data frame with model, sample, response, predictor importance and dependency 
  if (length(strsplit(ID, split = "_")[[1]]) < 3){ #samples not containing an underscore in the name
    RF_Imp.mt <- rbind(RF_Imp.mt,
                       c(ID,
                         strsplit(ID, split = "_")[[1]][1],
                         "",
                         strsplit(ID, split = "_")[[1]][2],
                         t(as.matrix(imp_plot$data$importance)),
                         cor(RF[[ID]][["observations"]][,"H1_TME"],RF[[ID]][["shaps"]][,"H1_TME"]),
                         cor(RF[[ID]][["observations"]][,"H13_TME"],RF[[ID]][["shaps"]][,"H13_TME"]),
                         cor(RF[[ID]][["observations"]][,"H3_TME"],RF[[ID]][["shaps"]][,"H3_TME"]),
                         cor(RF[[ID]][["observations"]][,"H5_TME"],RF[[ID]][["shaps"]][,"H5_TME"]),
                         cor(RF[[ID]][["observations"]][,"H6_TME"],RF[[ID]][["shaps"]][,"H6_TME"]),
                         cor(RF[[ID]][["observations"]][,"H7_TME"],RF[[ID]][["shaps"]][,"H7_TME"]),
                         cor(RF[[ID]][["observations"]][,"H9_TME"],RF[[ID]][["shaps"]][,"H9_TME"])
                         
                       ))
    
  } else { 
    RF_Imp.mt <- rbind(RF_Imp.mt,  #samples containing 1 underscore in the name
                       c(ID,
                         paste0(strsplit(ID, split = "_")[[1]][1],"_",strsplit(ID, split = "_")[[1]][2]) ,
                         "",
                         strsplit(ID, split = "_")[[1]][3],
                         t(as.matrix(imp_plot$data$importance)),
                         cor(RF[[ID]][["observations"]][,"H1_TME"],RF[[ID]][["shaps"]][,"H1_TME"]),
                         cor(RF[[ID]][["observations"]][,"H13_TME"],RF[[ID]][["shaps"]][,"H13_TME"]),
                         cor(RF[[ID]][["observations"]][,"H3_TME"],RF[[ID]][["shaps"]][,"H3_TME"]),
                         cor(RF[[ID]][["observations"]][,"H5_TME"],RF[[ID]][["shaps"]][,"H5_TME"]),
                         cor(RF[[ID]][["observations"]][,"H6_TME"],RF[[ID]][["shaps"]][,"H6_TME"]),
                         cor(RF[[ID]][["observations"]][,"H7_TME"],RF[[ID]][["shaps"]][,"H7_TME"]),
                         cor(RF[[ID]][["observations"]][,"H9_TME"],RF[[ID]][["shaps"]][,"H9_TME"])))
    
  }
  
}
#set column names to the importances
colnames(RF_Imp.mt) <- c("ID", "sample", "type", "response", "H1", "H13", "H3", "H5", "H6", "H7", "H9", "H1_cor", "H13_cor", "H3_cor", "H5_cor", "H6_cor", "H7_cor", "H9_cor")
#set columns as numeric
RF_Imp.mt[,c(5:18)] <- apply(RF_Imp.mt[,c(5:18)], 2, as.numeric)

# Add tumor type info
source("../utils/SamplesMetadata.R")
for (n in unique(RF_Imp.mt$sample)) {
  RF_Imp.mt$type[RF_Imp.mt$sample==n] <- annotation_site[n]
}
#set tumor type information as string/character
RF_Imp.mt$type <- as.character(RF_Imp.mt$type)
#only keep tumor type name
RF_Imp.mt$type <- sub(" .*", "", RF_Imp.mt$type)
#save importance matrix
write.table(RF_Imp.mt, file = "", sep = "\t")


### Extract R-squared
RF_R.mt <- data.frame(matrix(nrow = 0, ncol = 5))
#extract R-squared of each RF model
for (ID in names(RF)){
  #create data frame with model, sample, response, R-squared of the model 
  if (length(strsplit(ID, split = "_")[[1]]) < 3){ #samples not containing an underscore in the name
    RF_R.mt <- rbind(RF_R.mt,
                     c(ID,
                       strsplit(ID, split = "_")[[1]][1],
                       "",
                       strsplit(ID, split = "_")[[1]][2],
                       RF[[ID]][["r.squared"]]))
    
  } else { #samples containing 1 underscore in the name
    RF_R.mt <- rbind(RF_R.mt, 
                     c(ID,
                       paste0(strsplit(ID, split = "_")[[1]][1],"_",strsplit(ID, split = "_")[[1]][2]) ,
                       "",
                       strsplit(ID, split = "_")[[1]][3],
                       RF[[ID]][["r.squared"]]))
    
  }
  
}
#set column names to the R-squared
colnames(RF_R.mt) <- c("ID", "sample", "type", "response", "Rsquared")

RF_R.mt$Rsquared<- as.numeric(RF_R.mt$Rsquared)

# Add tumor type info
for (n in unique(RF_R.mt$sample)) {
  RF_R.mt$type[RF_R.mt$sample==n] <- annotation_site[n]
}
#set tumor type information as string/character
RF_R.mt$type <- as.character(RF_R.mt$type)
#only keep tumor type name
RF_R.mt$type <- sub(" .*", "", RF_R.mt$type)



#Import RF importances 
RF_Imp.mt <- read.delim("")
#Import RF R-squared
RF_R.mt <- read.delim("")
#Import data frame containing SCD metadata
df_meta <- read.delim("")
df_meta <- tibble::rownames_to_column(df_meta, "sample")

#Merge SCD values and Rsquared to the importance matrix
RF_Imp.mt <- merge(df_meta[,c("sample", "SCD")], RF_Imp.mt, by = "sample")
RF_Imp.mt <- merge(RF_R.mt[,c("ID", "Rsquared")], RF_Imp.mt, by = "ID")
#transpose the matrix and set rownames the first column
RF_Imp.mt.t <- data.frame(t(data.frame(RF_Imp.mt[,c(1,7:13)], row.names = 1)))
#compute the importance fraction for the heatmap
RF_Imp.mt.hmap <- sapply(RF_Imp.mt.t, function(x) {
  x = x / sum(x)
})
#set rownames
rownames(RF_Imp.mt.hmap) <- rownames(RF_Imp.mt.t)
#transpose again and convert to data frame
RF_Imp.mt.hmap <- data.frame(t(RF_Imp.mt.hmap))
#set rownames as a new column (ID)
RF_Imp.mt.hmap <- tibble::rownames_to_column(RF_Imp.mt.hmap, "ID")
#merge importance fraction and metadata and shapley dependencies
RF_Imp.mt.hmap <- merge(RF_Imp.mt.hmap, RF_Imp.mt[c(1:6, 14:20)], by= "ID")
RF_Imp.mt.hmap <- arrange(RF_Imp.mt.hmap, response,type)
#save tumor types names that will split the heatmap
split <- RF_Imp.mt.hmap$type
split <- factor(split, levels = unique(RF_Imp.mt.hmap$type))

#create a matrix for each information layer of the circos
RF_Imp <- data.matrix(data.frame(RF_Imp.mt.hmap[,c(1:8)], row.names = 1))
RF_Cor <- data.matrix(data.frame(RF_Imp.mt.hmap[,c(1,14:20)], row.names = 1))
RF_SCD <- data.matrix(data.frame(RF_Imp.mt.hmap[,c(1,11)], row.names = 1))
RF_response <- data.frame(RF_Imp.mt.hmap[,c(1,13)], row.names = 1)
RF_R <- data.frame(RF_Imp.mt.hmap[,c(1,9)], row.names = 1)
RF_sample <- data.frame(RF_Imp.mt.hmap[,c(1,10,12,13)], row.names = 1)
RF_sample <- arrange(RF_sample, type, response, sample) 

#number the samples within each tumor type
RF_sample$sample_ID <- ""
for (t in unique(RF_sample$type)){
  for (r in unique(RF_sample$response)){
    RF_sample$sample_ID[RF_sample$type ==t & RF_sample$response == r] <- 1:nrow(RF_sample[RF_sample$type ==t & RF_sample$response == r,])
  }
}
#arrange by response and tumor type
RF_sample <- arrange(RF_sample, response, type)
#remove the 3 first columns
RF_sample[1:3] <- NULL
colnames(RF_sample) <- "sample"

## Code to plot the detailed circos 
circos.clear()
circos.par(start.degree = 53, gap.after = c(3,3,3,3,3,3,3,3,30,2))

col_fun1 = colorRamp2(c(0, 0.5, 0.82), rev(paletteer_c("grDevices::Reds 3", 3)))
circos.heatmap(RF_Imp, split = split, col = col_fun1, track.height = 0.15,
               bg.border = "darkgreen", bg.lwd = 2, bg.lty = 3, show.sector.labels = TRUE,
               dend.side = "outside", clustering.method = "average", distance.method = "canberra")

circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 9) { # the last sector
    # predictors left
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -1.8,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -0.7,
                col = "#15CE59", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -3,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -1.8,
                col = "#95641A", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -4.2,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -3,
                col = "#CB3BBD", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -5.4,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -4.2,
                col = "#E6880D", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -6.6,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -5.4,
                col = "#000000", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -7.8,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -6.6,
                col = "#EE0F16", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -9,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -7.8,
                col = "#8E909B", border = NA)
    # importance annotation
    circos.rect(CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -0.5,
                CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -9,
                col = "blanchedalmond", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -2.5,
                "Fraction", cex = 0.4, facing = "bending")
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -5,
                "Feature", cex = 0.4, facing = "bending")
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -7.5,
                "Importances", cex = 0.38, facing = "bending")
    
    # predictors right
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -1.8,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -0.7,
                col = "#15CE59", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -3,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -1.8,
                col = "#95641A", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -4.2,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -3,
                col = "#CB3BBD", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -5.4,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -4.2,
                col = "#E6880D", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -6.6,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -5.4,
                col = "#000000", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -7.8,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -6.6,
                col = "#EE0F16", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -9,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -7.8,
                col = "#8E909B", border = NA)
    # direction
    # predictors left
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -11,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -9.8,
                col = "#15CE59", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -12.2,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -11,
                col = "#95641A", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -13.4,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -12.2,
                col = "#CB3BBD", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -14.5,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -13.3,
                col = "#E6880D", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -15.6,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -14.4,
                col = "#000000", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -16.7,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -15.5,
                col = "#EE0F16", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -17.8,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -16.6,
                col = "#8E909B", border = NA)
    
    # direction annotation
    circos.rect(CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -9.6,
                CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -17.8,
                col = "blanchedalmond", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -11,
                "Feature", cex = 0.39, facing = "bending")
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -13.5,
                "Dependency", cex = 0.35, facing = "bending")
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -16,
                "Direction", cex = 0.39, facing = "bending")
    
    # predictions right
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -11,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -9.8,
                col = "#15CE59", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -12.2,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -11,
                col = "#95641A", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -13.4,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -12.2,
                col = "#CB3BBD", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -14.5,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -13.3,
                col = "#E6880D", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -15.6,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -14.4,
                col = "#000000", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -16.7,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -15.5,
                col = "#EE0F16", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -17.8,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -16.6,
                col = "#8E909B", border = NA)
    
    # R2 and the rest
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -20.5,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -18.7,
                col = "blanchedalmond", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -19.2,
                "sample", cex = 0.5, facing = "bending")
    
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -23.1,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -21.3,
                col = "blanchedalmond", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -22,
                "R2", cex = 0.5, facing = "bending")
    
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -25.7,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -23.9,
                col = "blanchedalmond", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -24.5,
                "SCD", cex = 0.5, facing = "bending")
    
    
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -29.4,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -26.5,
                col = "blanchedalmond", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -28,
                "response", cex = 0.45, facing = "bending")
    
  }
}, bg.border = NA)


col_fun2 = colorRamp2(c(-1, 0, 1), c("blue", "white", "darkgreen"))
circos.heatmap(RF_Cor, col = col_fun2, track.height = 0.15)

col_sample <- structure(Seurat::DiscretePalette(15, palette = "polychrome") [c(1:3, 5, 6:10)], names = unique(RF_sample$sample))
circos.heatmap(RF_sample, col = col_sample, track.height = 0.03)

col_fun3 = colorRamp2(c(0, 1), c("white", "black"))
circos.heatmap(RF_R, col = col_fun3, track.height = 0.03)

col_fun4 = colorRamp2(c(0.8, 1), c("white", "black"))
circos.heatmap(RF_SCD, col = col_fun4, track.height = 0.03)

col_response <- structure(c("#71189E", "#05F3EB", "#890269", "#701717", "#4969D4", "#132892"), names = unique(RF_response$response))
circos.heatmap(RF_response, col = col_response, track.height = 0.05)

### legends
lgd_Imp = Legend(title = "Fraction of Importances", col_fun = col_fun1)
lgd_dir = Legend(title = "Feature Dependency", col_fun =  col_fun2)
lgd_R = Legend(title = "R2", col_fun = col_fun3)
lgd_SCD = Legend(title = "Spatial Continuity Degree", col_fun = col_fun4)
lgd_resp = Legend(title = "Response", at = names(col_response),
                  legend_gp = gpar(fill = col_response))
lgd_sample = Legend(title = "Sample", at = names(col_sample),
                    legend_gp = gpar(fill = col_sample))

library(gridBase)
circos.clear()
dev.off()
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region

pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)

upViewport()

h = dev.size()[1]
lgd_list = packLegend(lgd_Imp, lgd_dir, lgd_R, lgd_SCD,
                      lgd_resp,lgd_sample, max_height = unit(0.9*h, "inch"))
draw(lgd_list, x = circle_size, just = "left")

#Create the circos links
RF_Imp_direction <- RF_Imp.mt
RF_Imp_direction[,c(2,4)] <- NULL
RF_Imp_direction$predictor <- ""
RF_Imp_direction$direction <- ""

for (id in unique(RF_Imp_direction$ID)){
  top_pred <- names(which.max(RF_Imp_direction[RF_Imp_direction$ID == id,5:11]))
  if (RF_Imp_direction[RF_Imp_direction$ID == id, paste0(top_pred, "_cor")] >= 0.6) {
    RF_Imp_direction[RF_Imp_direction$ID == id,]$predictor <- top_pred
    RF_Imp_direction[RF_Imp_direction$ID == id,]$direction <- "pos"
  } else if (RF_Imp_direction[RF_Imp_direction$ID == id, paste0(top_pred, "_cor")] <= - 0.6){
    RF_Imp_direction[RF_Imp_direction$ID == id,]$predictor <- top_pred
    RF_Imp_direction[RF_Imp_direction$ID == id,]$direction <- "neg"
  } else {
    RF_Imp_direction[RF_Imp_direction$ID == id,]$predictor <- top_pred
    RF_Imp_direction[RF_Imp_direction$ID == id,]$direction <- "nonlinear"
  }
}
table(RF_Imp_direction$direction)
table(RF_Imp_direction$direction, RF_Imp_direction$predictor)

write.table(RF_Imp_direction, "")

RF_Imp_direction <- filter(RF_Imp_direction, direction != "nonlinear")
RF_Imp_direction$link <- ""

for (s in unique(RF_Imp_direction$sample)) {
  for (p in paste0("H", c(1,3,5,6,7,9,13))){
    print(p)
    if (p %in% RF_Imp_direction[RF_Imp_direction$sample == s,]$predictor) {
      if ((table(RF_Imp_direction[RF_Imp_direction$sample == s,]$predictor == p)["TRUE"] >= 2) &
          (length(unique(RF_Imp_direction[RF_Imp_direction$sample == s & RF_Imp_direction$predictor == p,]$direction)) > 1)) {
        RF_Imp_direction[RF_Imp_direction$sample == s & RF_Imp_direction$predictor == p,]$link <- paste0(p, "_link")
      }
    }
  }
}

RF_Imp_direction <- merge(tibble::rownames_to_column(data.frame(RF_Imp), "ID"), RF_Imp_direction[,c("ID", "sample", "link")], by = "ID", all=T)
unique(RF_Imp_direction$link)

# H3 only positive
id <- RF_Imp_direction$ID[RF_Imp_direction$predictor == "H3" & RF_Imp_direction$direction == "pos"]
idx <- which(rownames(RF_Imp) %in% id)
tmp <- t(combn(idx, 2))
df_link_H3 <- data.frame(tmp)
colnames(df_link_H3) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H3))) {
  circos.heatmap.link(df_link_H3$from_index[i],
                      df_link_H3$to_index[i],
                      col = c("#CB3BBD"), lwd = 0.5)
}


# H1 only positive
id <- RF_Imp_direction$ID[RF_Imp_direction$predictor == "H1" & RF_Imp_direction$direction == "pos"]
idx <- which(rownames(RF_Imp) %in% id)
tmp <- t(combn(idx, 2))
df_link_H1 <- data.frame(tmp)
colnames(df_link_H1) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H1))) {
  circos.heatmap.link(df_link_H1$from_index[i],
                      df_link_H1$to_index[i],
                      col = c("#15CE59"), lwd = 0.5)
}

# H7 only positive
id <- RF_Imp_direction$ID[RF_Imp_direction$predictor == "H7" & RF_Imp_direction$direction == "pos"]
idx <- which(rownames(RF_Imp) %in% id)
tmp <- t(combn(idx, 2))
df_link_H7 <- data.frame(tmp)
colnames(df_link_H7) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H7))) {
  circos.heatmap.link(df_link_H7$from_index[i],
                      df_link_H7$to_index[i],
                      col = c("#EE0F16"), lwd = 0.5)
}


# H13 only positive
id <- RF_Imp_direction$ID[RF_Imp_direction$predictor == "H13" & RF_Imp_direction$direction == "pos"]
idx <- which(rownames(RF_Imp) %in% id)
tmp <- t(combn(idx, 2))
df_link_H13 <- data.frame(tmp)
colnames(df_link_H13) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H13))) {
  circos.heatmap.link(df_link_H13$from_index[i],
                      df_link_H13$to_index[i],
                      col = c("#95641A"), lwd = 0.5)
}


# H9 only positive
id <- RF_Imp_direction$ID[RF_Imp_direction$predictor == "H9" & RF_Imp_direction$direction == "pos"]
idx <- which(rownames(RF_Imp) %in% id)
tmp <- t(combn(idx, 2))
df_link_H9 <- data.frame(tmp)
colnames(df_link_H9) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H9))) {
  circos.heatmap.link(df_link_H9$from_index[i],
                      df_link_H9$to_index[i],
                      col = c("#8E909B"), lwd = 0.5)
}


# H5 only positive
id <- RF_Imp_direction$ID[RF_Imp_direction$predictor == "H5" & RF_Imp_direction$direction == "pos"]
idx <- which(rownames(RF_Imp) %in% id)
tmp <- t(combn(idx, 2))
df_link_H5 <- data.frame(tmp)
colnames(df_link_H5) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H5))) {
  circos.heatmap.link(df_link_H5$from_index[i],
                      df_link_H5$to_index[i],
                      col = c("#E6880D"), lwd = 0.5)
}

# H6 only positive
id <- RF_Imp_direction$ID[RF_Imp_direction$predictor == "H6" & RF_Imp_direction$direction == "pos"]
idx <- which(rownames(RF_Imp) %in% id)
tmp <- t(combn(idx, 2))
df_link_H6 <- data.frame(tmp)
colnames(df_link_H6) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H6))) {
  circos.heatmap.link(df_link_H6$from_index[i],
                      df_link_H6$to_index[i],
                      col = c("#000000"), lwd = 0.5)
}

write.table(RF_Imp_direction, "")


################# General circos with boxplots ###################

# get the approporitate importance fractions dataframe
library(reshape2)
RF_Imp.t <- RF_Imp.mt.hmap[,c(1:8, 10, 12, 13)]
RF_Imp.t$response <- factor(RF_Imp.t$response, levels = c("H0", "H11", "H12", "H2", "H4", "H8"))
RF_Imp.t <- arrange(RF_Imp.t, response, type)
RF_Imp.t <- melt(RF_Imp.t[,c(1:8)], id.vars = "ID")
RF_Imp.t <- merge(RF_Imp.t, RF_Imp.mt.hmap[,c("ID", "type", "sample", "response")], by= "ID")
RF_Imp.t$ID <- paste0(RF_Imp.t$sample, "_", RF_Imp.t$variable)
RF_Imp.t$sample <- NULL

# wide
library(tidyr)
RF_Imp.t.wide <- spread(RF_Imp.t, response, value)

RF_Imp.t.wide <- arrange(RF_Imp.t.wide, variable, type)

RF_Imp.t.wide <- RF_Imp.t.wide[,-1] %>% group_by(type, variable) %>% summarise_all(mean)

# get the appropriate importance directions dataframe
library(reshape2)
RF_Dir.t <- RF_Imp.mt.hmap[,c(1, 14:20, 10, 12, 13)]
RF_Dir.t$response <- factor(RF_Dir.t$response, levels = c("H0", "H11", "H12", "H2", "H4", "H8"))
RF_Dir.t <- arrange(RF_Dir.t, response, type)
RF_Dir.t <- melt(RF_Dir.t[,c(1:8)], id.vars = "ID")
RF_Dir.t <- merge(RF_Dir.t, RF_Imp.mt.hmap[,c("ID", "type", "sample", "response")], by= "ID")

RF_Dir.t$ID <- paste0(RF_Dir.t$sample, "_", RF_Dir.t$variable)
RF_Dir.t$sample <- NULL

# wide
library(tidyr)
RF_Dir.t.wide <- spread(RF_Dir.t, response, value)

RF_Dir.t.wide <- arrange(RF_Dir.t.wide, variable, type)

RF_Dir.t.wide <- RF_Dir.t.wide[,-1] %>% group_by(type, variable) %>% summarise_all(median)

## plot circos with boxplots
split <- RF_Dir.t.wide$type
split <- factor(split, levels = unique(RF_Dir.t.wide$type))

RF_Imp.t <- data.matrix(RF_Imp.t.wide[,c(3:8)])
RF_Dir.t <- data.matrix(RF_Dir.t.wide[,c(3:8)])                        

circos.clear()

circos.par(start.degree = 42, gap.after = c(3,3,3,3,3,3,3,3,30,2))

col_fun2 = colorRamp2(c(-1, 0, 1), c("blue", "white", "darkgreen"))
circos.heatmap(RF_Dir.t, split = split, col = col_fun2, track.height = 0.2,
               bg.border = "darkred", bg.lwd = 2, bg.lty = 3, show.sector.labels = TRUE, cluster=T,
               dend.side = "outside", clustering.method = "average", distance.method = "canberra")

circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 9) { # the last sector
    # predictors left
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -2.9,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -1,
                col = "#71189E", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -4.8,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -2.9,
                col = "#05F3EB", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -6.7,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -4.8,
                col = "#890269", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -8.6,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -6.7,
                col = "#701717", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -10.5,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -8.6,
                col = "#4969D4", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -12.5,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -10.4,
                col = "#132892", border = NA)
    # importance annotation
    circos.rect(CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -12.5,
                CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -1,
                col = "blanchedalmond", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -2,
                "Cancer", cex = 0.4, facing = "bending")
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -5,
                "Response", cex = 0.4, facing = "bending")
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -8,
                "Dependency", cex = 0.38, facing = "bending")
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -11,
                "Direction", cex = 0.38, facing = "bending")
    
    # predictors right
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -2.9,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -1,
                col = "#71189E", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -4.8,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -2.9,
                col = "#05F3EB", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -6.7,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -4.8,
                col = "#890269", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -8.6,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -6.7,
                col = "#701717", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -10.5,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -8.6,
                col = "#4969D4", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -12.5,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -10.4,
                col = "#132892", border = NA)
    
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -15.6,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -13.4,
                col = "blanchedalmond", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -14,
                "predictor", cex = 0.5, facing = "bending")
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -28,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -16.3,
                col = "blanchedalmond", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -17.3,
                "Mean", cex = 0.5, facing = "bending")
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -20,
                "Fraction", cex = 0.5, facing = "bending")
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -23,
                "Feature", cex = 0.5, facing = "bending")
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -26,
                "Importance", cex = 0.45, facing = "bending")
  }
}, bg.border = NA)


col_pred <- structure(c("#15CE59", "#95641A", "#CB3BBD", "#E6880D", "#000000", "#EE0F16", "#8E909B"), names = unique(as.character(RF_Dir.t.wide$variable)))
circos.heatmap(as.matrix(RF_Dir.t.wide[,"variable"]), col = col_pred, track.height = 0.03)

circos.track(ylim = range(RF_Imp.t), panel.fun = function(x, y) {
  m = RF_Imp.t[CELL_META$subset, 1:5, drop = FALSE]
  m = m[CELL_META$row_order, , drop = FALSE]
  n = nrow(m)
  # circos.boxplot is applied on matrix columns, so here we transpose it.
  circos.boxplot(t(m), pos = 1:n - 0.5, pch = 16, cex = 0.3)
  circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey")
}, cell.padding = c(0.02, 0, 0.02, 0))

### legends
library(ComplexHeatmap)
lgd_dir = Legend(title = "Feature Dependency", col_fun =  col_fun2)

library(gridBase)
circos.clear()
dev.off()
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region

pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)

upViewport()

h = dev.size()[1]
lgd_list = packLegend(lgd_dir, max_height = unit(0.9*h, "inch"))
draw(lgd_list, x = circle_size, just = "left")

### Add links for importances
avg_pred <- RF_Imp.mt.hmap[,c(2:8, 12)]
avg_pred <- avg_pred %>% group_by(type) %>% summarise_all(mean)
avg_pred$top_pred <- ""

for(r in 1:nrow(avg_pred)){
  avg_pred[r,]$top_pred<- names(which.max(avg_pred[r,-1]))
}

avg_pred$variable <- paste0(avg_pred$top_pred, "_cor")

RF_Dir.t.wide <- merge(RF_Dir.t.wide, avg_pred[,c("type", "top_pred", "variable")], by = c("type", "variable"), all = T)
colnames(RF_Dir.t.wide)[9] <- "link"

### Add links for direction
RF_Dir.t.wide$link <- ""

for (t in unique(RF_Dir.t.wide$type)){
  for (p in unique(RF_Dir.t.wide$variable)){
    if (rowSums(filter(RF_Dir.t.wide, type == t & variable == p)[,c(3:8)] >= 0.2) >= 4){
      RF_Dir.t.wide[RF_Dir.t.wide$type == t & RF_Dir.t.wide$variable == p,]$link <- paste0(p, "_link")
    }
  }
}

table(RF_Dir.t.wide$link)

# H9 link
id <- rownames(RF_Dir.t.wide)[RF_Dir.t.wide$link == "H9_cor_link" & !is.na(RF_Dir.t.wide$link)]
idx <- which(rownames(as.data.frame(RF_Dir.t)) %in% id)
tmp <- t(combn(idx, 2))
df_link_H9 <- data.frame(tmp)
colnames(df_link_H9) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H9))) {
  circos.heatmap.link(df_link_H9$from_index[i],
                      df_link_H9$to_index[i],
                      col = c("#8E909B"), lwd= 3)
}

# H1 link
id <- rownames(RF_Dir.t.wide)[RF_Dir.t.wide$link == "H1_cor_link" & !is.na(RF_Dir.t.wide$link)]
idx <- which(rownames(as.data.frame(RF_Dir.t)) %in% id)
tmp <- t(combn(idx, 2))
df_link_H1 <- data.frame(tmp)
colnames(df_link_H1) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H1))) {
  circos.heatmap.link(df_link_H1$from_index[i],
                      df_link_H1$to_index[i],
                      col = c("#15CE59"), lwd= 3)
}


# H3 link
id <- rownames(RF_Dir.t.wide)[RF_Dir.t.wide$link == "H3_cor_link" & !is.na(RF_Dir.t.wide$link)]
idx <- which(rownames(as.data.frame(RF_Dir.t)) %in% id)
tmp <- t(combn(idx, 2))
df_link_H3 <- data.frame(tmp)
colnames(df_link_H3) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H3))) {
  circos.heatmap.link(df_link_H3$from_index[i],
                      df_link_H3$to_index[i],
                      col = c("#CB3BBD"), lwd= 3)
}

# H6 link
id <- rownames(RF_Dir.t.wide)[RF_Dir.t.wide$link == "H6_cor_link" & !is.na(RF_Dir.t.wide$link)]
idx <- which(rownames(as.data.frame(RF_Dir.t)) %in% id)
tmp <- t(combn(idx, 2))
df_link_H6 <- data.frame(tmp)
colnames(df_link_H6) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H6))) {
  circos.heatmap.link(df_link_H6$from_index[i],
                      df_link_H6$to_index[i],
                      col = c("#000000"), lwd= 3)
}

# H5 link
id <- rownames(RF_Dir.t.wide)[RF_Dir.t.wide$link == "H5_cor_link" & !is.na(RF_Dir.t.wide$link)]
idx <- which(rownames(as.data.frame(RF_Dir.t)) %in% id)
tmp <- t(combn(idx, 2))
df_link_H5 <- data.frame(tmp)
colnames(df_link_H5) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H5))) {
  circos.heatmap.link(df_link_H5$from_index[i],
                      df_link_H5$to_index[i],
                      col = c("#E6880D"), lwd= 3)
}


# H13 link
id <- rownames(RF_Dir.t.wide)[RF_Dir.t.wide$link == "H13_cor_link" & !is.na(RF_Dir.t.wide$link)]
idx <- which(rownames(as.data.frame(RF_Dir.t)) %in% id)
tmp <- t(combn(idx, 2))
df_link_H13 <- data.frame(tmp)
colnames(df_link_H13) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H13))) {
  circos.heatmap.link(df_link_H13$from_index[i],
                      df_link_H13$to_index[i],
                      col = c("#95641A"), lwd= 3)
}

# H7 link
id <- rownames(RF_Dir.t.wide)[RF_Dir.t.wide$link == "H7_cor_link" & !is.na(RF_Dir.t.wide$link)]
idx <- which(rownames(as.data.frame(RF_Dir.t)) %in% id)
tmp <- t(combn(idx, 2))
df_link_H7 <- data.frame(tmp)
colnames(df_link_H7) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H7))) {
  circos.heatmap.link(df_link_H7$from_index[i],
                      df_link_H7$to_index[i],
                      col = c("#EE0F16"),lwd= 3)
}



#to be inspected later#

##############
RF_Imp_direction <- merge(RF_Imp_direction, df_meta, by = "sample", all =T)
TME_link <- RF_Imp_direction[,c("link", "Buffer", "Cancer", "TME")]
TME_link <- melt(TME_link, id.vars = "link")

ggboxplot(filter(TME_link, link != "" & is.na(link) == F), x= "link", y = "value", facet.by = "variable") + theme(axis.text.x = element_text(angle = -45, hjust = 0))


# link and sample (or type)
TME_link <- RF_Imp_direction[,c("link", "sample", "direction", "response")]
TME_link <- filter(TME_link, direction != "nonlinear" & link != "")
#TME_link[,c("direction", "response")] <- NULL
TME_link <- TME_link[duplicated(TME_link)==F,]

TME_link <- TME_link[,c("link", "sample")]
TME_link <- TME_link[duplicated(TME_link)==F,]

rowSums(table(TME_link$link, TME_link$sample))
unique(TME_link$sample)

x = list("H1_link" = filter(TME_link, link == "H1_link")$sample,
         "H3_link" = filter(TME_link, link == "H3_link")$sample,
         "H5_link" = filter(TME_link, link == "H5_link")$sample,
         "H6_link" = filter(TME_link, link == "H6_link")$sample,
         "H7_link" = filter(TME_link, link == "H7_link")$sample,
         "H9_link" = filter(TME_link, link == "H9_link")$sample,
         "H13_link" = filter(TME_link, link == "H13_link")$sample)

x = list("H1_link" = filter(TME_link, link == "H1_link")$type,
         "H3_link" = filter(TME_link, link == "H3_link")$type,
         "H5_link" = filter(TME_link, link == "H5_link")$type,
         "H6_link" = filter(TME_link, link == "H6_link")$type,
         "H7_link" = filter(TME_link, link == "H7_link")$type,
         "H9_link" = filter(TME_link, link == "H9_link")$type,
         "H13_link" = filter(TME_link, link == "H13_link")$type
)

#ggVennDiagram(x, label = "count", edge_size = 2, label_geom = "text") + scale_fill_distiller(palette = "RdYlBu")
library(UpSetR)
upset(fromList(x), nsets = 7, sets.x.label = "Total = 26 samples", text.scale = 2)


###
TME_link <- RF_Imp_direction[,c("link", "sample", "direction", "response")]
TME_link <- filter(TME_link, direction != "nonlinear" & link != "")
TME_link <- TME_link[duplicated(TME_link)==F,]


sum(rowSums(table(TME_link$link, TME_link$sample)))
colSums(table(TME_link$link, TME_link$type))

# samle and response (by link)

# H3
TME_link_H3 <- filter(TME_link, link == "H3_link")
TME_link_H3.tbl <- data.frame(table(TME_link_H3$sample, TME_link_H3$response))
colnames(TME_link_H3.tbl)[1:2] <- c("sample", "response")
TME_link_H3.tbl <- merge(TME_link_H3.tbl, TME_link_H3[,-1], by = c("sample", "response"), all=T)
TME_link_H3.tbl$link <- "H3_link"
# H5
TME_link_H5 <- filter(TME_link, link == "H5_link")
TME_link_H5.tbl <- data.frame(table(TME_link_H5$sample, TME_link_H5$response))
colnames(TME_link_H5.tbl)[1:2] <- c("sample", "response")
TME_link_H5.tbl <- merge(TME_link_H5.tbl, TME_link_H5[,-1], by = c("sample", "response"), all=T)
TME_link_H5.tbl$link <- "H5_link"
# H6
TME_link_H6 <- filter(TME_link, link == "H6_link")
TME_link_H6.tbl <- data.frame(table(TME_link_H6$sample, TME_link_H6$response))
colnames(TME_link_H6.tbl)[1:2] <- c("sample", "response")
TME_link_H6.tbl <- merge(TME_link_H6.tbl, TME_link_H6[,-1], by = c("sample", "response"), all=T)
TME_link_H6.tbl$link <- "H6_link"
# H7
TME_link_H7 <- filter(TME_link, link == "H7_link")
TME_link_H7.tbl <- data.frame(table(TME_link_H7$sample, TME_link_H7$response))
colnames(TME_link_H7.tbl)[1:2] <- c("sample", "response")
TME_link_H7.tbl <- merge(TME_link_H7.tbl, TME_link_H7[,-1], by = c("sample", "response"), all=T)
TME_link_H7.tbl$link <- "H7_link"
# H9
TME_link_H9 <- filter(TME_link, link == "H9_link")
TME_link_H9.tbl <- data.frame(table(TME_link_H9$sample, TME_link_H9$response))
colnames(TME_link_H9.tbl)[1:2] <- c("sample", "response")
TME_link_H9.tbl <- merge(TME_link_H9.tbl, TME_link_H9[,-1], by = c("sample", "response"), all=T)
TME_link_H9.tbl$link <- "H9_link"
# H13
TME_link_H13 <- filter(TME_link, link == "H13_link")
TME_link_H13.tbl <- data.frame(table(TME_link_H13$sample, TME_link_H13$response))
colnames(TME_link_H13.tbl)[1:2] <- c("sample", "response")
TME_link_H13.tbl <- merge(TME_link_H13.tbl, TME_link_H13[,-1], by = c("sample", "response"), all=T)
TME_link_H13.tbl$link <- "H13_link"
# H1
TME_link_H1 <- filter(TME_link, link == "H1_link")
TME_link_H1.tbl <- data.frame(table(TME_link_H1$sample, TME_link_H1$response))
colnames(TME_link_H1.tbl)[1:2] <- c("sample", "response")
TME_link_H1.tbl <- merge(TME_link_H1.tbl, TME_link_H1[,-1], by = c("sample", "response"), all=T)
TME_link_H1.tbl$link <- "H1_link"

TME_link.all <- rbind(TME_link_H1.tbl, TME_link_H3.tbl, TME_link_H5.tbl, TME_link_H6.tbl, TME_link_H7.tbl, TME_link_H9.tbl, TME_link_H13.tbl)


col_direction <- structure(c("darkblue", "darkgreen"), names = c("neg", "pos"))
ggbarplot(TME_link.all, x= "response", y = "Freq", fill = "direction", palette = col_direction, ylab = "nSamples") + facet_wrap(~link, nrow = 5, ncol = 7) + theme(axis.text.x = element_text(angle = -45, hjust = 0))


TME_link_H10 <- filter(TME_link.all, response == "H10" & Freq > 0)
TME_link_H10 <- merge(TME_link_H10, filter(HallmarkActivity_labels_Cancer, variable == "H10")[,c("sample", "diff", "label")], by = "sample")

# check how H10 is up or down in all samples where it is in the links
table(TME_link_H10$label)
table(TME_link_H10$label, TME_link_H10$link)
table(TME_link_H10$label, TME_link_H10$direction)


# get the intersect of samples between the two links: cancer and TME
unique(TME_link$sample)
unique(cancer_link$sample)

common_samples <- unique(cancer_link$sample)[unique(cancer_link$sample) %in% unique(TME_link$sample) == T]

############### boxplots #################
###

RF_Imp_direction <- RF_Imp.mt
RF_Imp_direction[,c(2,4)] <- NULL
RF_Imp_direction$predictor <- ""
RF_Imp_direction$direction <- ""

for (id in unique(RF_Imp_direction$ID)){
  top_pred <- names(which.max(RF_Imp_direction[RF_Imp_direction$ID == id,5:11]))
  if (RF_Imp_direction[RF_Imp_direction$ID == id, paste0(top_pred, "_cor")] >= 0.3) {
    RF_Imp_direction[RF_Imp_direction$ID == id,]$predictor <- top_pred
    RF_Imp_direction[RF_Imp_direction$ID == id,]$direction <- "pos"
  } else if (RF_Imp_direction[RF_Imp_direction$ID == id, paste0(top_pred, "_cor")] <= - 0.3){
    RF_Imp_direction[RF_Imp_direction$ID == id,]$predictor <- top_pred
    RF_Imp_direction[RF_Imp_direction$ID == id,]$direction <- "neg"
  } else {
    RF_Imp_direction[RF_Imp_direction$ID == id,]$predictor <- top_pred
    RF_Imp_direction[RF_Imp_direction$ID == id,]$direction <- "nonlinear"
  }
}


freq_pan <- data.frame(table(RF_Imp_direction$predictor, RF_Imp_direction$response))
#ggboxplot(freq_pan, x="Var1", y = "Freq") + stat_compare_means(comparisons = list(c("H11", "H8"), c("H10", "H12"),
#                                                                                  c("H11", "H12"), c("H10", "H12")), method = "t.test")
color_codes <-  c(H1= "#15CE59",
                  H2 = "#701717",
                  H3 = "#CB3BBD",
                  H4 = "#4969D4",
                  H5 = "#E6880D",
                  H6 = "#000000",
                  H7 = "#EE0F16",
                  H8 = "#132892",
                  H9 = "#8E909B",
                  H10 = "#71189E",
                  H11 = "#05F3EB",
                  H12 = "#890269",
                  H13 = "#95641A")


freq_pan <- arrange(freq_pan, desc(Freq))
freq_pan$Var1 <- factor(freq_pan$Var1, levels = unique(freq_pan$Var1))

freq_pan$color <- freq_pan$Var1
freq_pan$color <- as.character(freq_pan$color)

for (hallmark in paste0("H", c(1,3,5,6,7,9,13))) {
  freq_pan$color[freq_pan$color == hallmark] <- color_codes[[hallmark]]
}

ggboxplot(freq_pan, x="Var1", y = "Freq", fill = "Var1", palette = unique(freq_pan$color)) + stat_compare_means(comparisons = list(c("H3", "H13")))              


## make a table of fractions distribution (take the median)
freq_pan_2 <- melt(RF_Imp.mt.hmap[,c(2:8,13)], id.vars = "response")
freq_pan_2$dir <- ""
for (pred in paste0("H", c(1,3,5,6,7,9,13))) {
  for (resp in paste0("H", c(10,11,12,2,4,8))){
    freq_pan_2$dir[freq_pan_2$variable == pred & freq_pan_2$response == resp] <- median(filter(RF_Imp_direction[,c("response", paste0(pred, "_cor"))], response == resp)[,paste0(pred, "_cor")])
  }
}
freq_pan_2$dir <- as.numeric(freq_pan_2$dir)
freq_pan_2$variable <- factor(freq_pan_2$variable, levels = c("H3", "H1", "H9", "H6", "H7", "H13", "H5"))
freq_pan_2 <- arrange(freq_pan_2, variable)

p <-ggplot(freq_pan_2, aes(x=variable, y=value, fill = dir, group = dir)) +
  geom_boxplot()+
  #coord_cartesian(ylim = c(0, 0.4)) +
  #geom_dotplot(binaxis='y', stackdir='center',binpositions="all", stackgroups=F, dotsize = 0.15)+
  scale_fill_gradient2(
    low = muted("blue"),
    mid = "white",
    high = muted("darkgreen"),
    midpoint = 0,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
  ) + facet_wrap(~response, nrow = 1)

p


## use for stats
ggboxplot(freq_pan_2, x="variable", y = "value") +
  stat_compare_means(comparisons = list(c("H2", "H8")),  
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                        symbols = c("****", "***", "**", "*", "ns")))


ggboxplot(freq_pan_2, x="variable", y = "value") +
  stat_compare_means(ref.group = c("H3"))

ggboxplot(freq_pan_2, x="variable", y = "value") +
  stat_compare_means()

## use this for plot
library(ggpubr)
ggplot(freq_pan_2, aes(x = variable, y = value, fill = dir)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch = 21, position = position_jitterdodge()) +
  #coord_cartesian(ylim = c(0, 0.68)) +
  scale_fill_gradient2(
    low = muted("blue"),
    mid = "white",
    high = muted("darkgreen"),
    midpoint = 0,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
  ) + theme_classic() + ylab("Feature Importance Fraction")  +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_blank(),
        legend.title = element_blank()
  ) +
  stat_compare_means(ref.group = c("H3"), size = 4.5, label = "p.format",family = "Times New Roman", label.y = 0.9)+
  stat_compare_means(family = "Times New Roman", size = 5, label.y = 1)
