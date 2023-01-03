##################################################
## Project: Cancer Hallmarks
## Script purpose: Plot results of Random Forest models for TME Hallmarks
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
############################## prediction of TME compartment ###################

RF <- list()
#grab each Random Forest shapley model and put it in the same list
for(file in list.files("")) {
  print(file)
  RF[[strsplit(file, split = ".rds")[[1]]]] <- readRDS("")
}

### Extract Importances
RF_Imp.mt <- data.frame(matrix(nrow = 0, ncol = 16))

#extract importances of each RF model 
for (ID in names(RF)){
  #barplot of importances
  imp_plot <- plot_feature_importance(RF[[ID]], max_vars = 6)
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
                         cor(RF[[ID]][["observations"]][,"H10_cancer"],RF[[ID]][["shaps"]][,"H10_cancer"]),
                         cor(RF[[ID]][["observations"]][,"H11_cancer"],RF[[ID]][["shaps"]][,"H11_cancer"]),
                         cor(RF[[ID]][["observations"]][,"H12_cancer"],RF[[ID]][["shaps"]][,"H12_cancer"]),
                         cor(RF[[ID]][["observations"]][,"H2_cancer"],RF[[ID]][["shaps"]][,"H2_cancer"]),
                         cor(RF[[ID]][["observations"]][,"H4_cancer"],RF[[ID]][["shaps"]][,"H4_cancer"]),
                         cor(RF[[ID]][["observations"]][,"H8_cancer"],RF[[ID]][["shaps"]][,"H8_cancer"])
                       ))
    
  } else { #samples containing 1 underscore in the name
    RF_Imp.mt <- rbind(RF_Imp.mt,
                       c(ID,
                         paste0(strsplit(ID, split = "_")[[1]][1],"_",strsplit(ID, split = "_")[[1]][2]) ,
                         "",
                         strsplit(ID, split = "_")[[1]][3],
                         t(as.matrix(imp_plot$data$importance)),
                         cor(RF[[ID]][["observations"]][,"H10_cancer"],RF[[ID]][["shaps"]][,"H10_cancer"]),
                         cor(RF[[ID]][["observations"]][,"H11_cancer"],RF[[ID]][["shaps"]][,"H11_cancer"]),
                         cor(RF[[ID]][["observations"]][,"H12_cancer"],RF[[ID]][["shaps"]][,"H12_cancer"]),
                         cor(RF[[ID]][["observations"]][,"H2_cancer"],RF[[ID]][["shaps"]][,"H2_cancer"]),
                         cor(RF[[ID]][["observations"]][,"H4_cancer"],RF[[ID]][["shaps"]][,"H4_cancer"]),
                         cor(RF[[ID]][["observations"]][,"H8_cancer"],RF[[ID]][["shaps"]][,"H8_cancer"])))
    
  }
  
}
#set column names to the importances
colnames(RF_Imp.mt) <- c("ID", "sample", "type", "response", "H10", "H11", "H12", "H2", "H4", "H8", "H10_cor", "H11_cor", "H12_cor", "H2_cor", "H4_cor", "H8_cor")
#set columns as numeric
RF_Imp.mt[,c(5:16)] <- apply(RF_Imp.mt[,c(5:16)], 2, as.numeric)

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
RF_Imp.mt.t <- data.frame(t(data.frame(RF_Imp.mt[,c(1,7:12)], row.names = 1)))
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
RF_Imp.mt.hmap <- merge(RF_Imp.mt.hmap, RF_Imp.mt[c(1:6, 13:18)], by= "ID")
RF_Imp.mt.hmap <- arrange(RF_Imp.mt.hmap, response,type)
#save tumor types names that will split the heatmap
split <- RF_Imp.mt.hmap$type
split <- factor(split, levels = unique(RF_Imp.mt.hmap$type))

#create a matrix for each information layer of the circos
RF_Imp <- data.matrix(data.frame(RF_Imp.mt.hmap[,c(1:7)], row.names = 1))
RF_Cor <- data.matrix(data.frame(RF_Imp.mt.hmap[,c(1,13:18)], row.names = 1))
RF_SCD <- data.matrix(data.frame(RF_Imp.mt.hmap[,c(1,10)], row.names = 1))
RF_response <- data.frame(RF_Imp.mt.hmap[,c(1,12)], row.names = 1)
RF_R <- data.frame(RF_Imp.mt.hmap[,c(1,8)], row.names = 1)
RF_sample <- data.frame(RF_Imp.mt.hmap[,c(1,9,11,12)], row.names = 1)
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
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -1.9,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -0.7,
                col = "#71189E", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -3,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -1.8,
                col = "#05F3EB", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -4.1,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -2.9,
                col = "#890269", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -5.2,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -4,
                col = "#701717", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -6.3,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -5.1,
                col = "#4969D4", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -7.4,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -6.2,
                col = "#132892", border = NA)
    
    # importance annotation
    circos.rect(CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -0.6,
                CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -7.5,
                col = "blanchedalmond", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -1.5,
                "Fraction", cex = 0.4, facing = "bending")
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -3.5,
                "Feature", cex = 0.4, facing = "bending")
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -5.5,
                "Importances", cex = 0.38, facing = "bending")
    
    # predictors right
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -1.9,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -0.7,
                col = "#71189E", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -3,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -1.8,
                col = "#05F3EB", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -4.1,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -2.9,
                col = "#890269", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -5.2,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -4,
                col = "#701717", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -6.3,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -5.1,
                col = "#4969D4", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -7.4,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -6.2,
                col = "#132892", border = NA)
    
    # direction
    # predictors left
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -9.4,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -8.2,
                col = "#71189E", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -10.5,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -9.3,
                col = "#05F3EB", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -11.6,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -10.4,
                col = "#890269", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -12.7,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -11.5,
                col = "#701717", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -13.8,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -12.6,
                col = "#4969D4", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -14.9,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -13.7,
                col = "#132892", border = NA)
    # direction annotation
    circos.rect(CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -8.1,
                CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -15,
                col = "blanchedalmond", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -9,
                "Feature", cex = 0.39, facing = "bending")
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -11,
                "Dependency", cex = 0.35, facing = "bending")
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -13,
                "Direction", cex = 0.39, facing = "bending")
    
    # predictors right
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -9.4,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -8.2,
                col = "#71189E", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -10.5,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -9.3,
                col = "#05F3EB", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -11.6,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -10.4,
                col = "#890269", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -12.7,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -11.5,
                col = "#701717", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -13.8,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -12.6,
                col = "#4969D4", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -14.9,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -13.7,
                col = "#132892", border = NA)
    
    # R2 and rest
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -17.2,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -15.7,
                col = "blanchedalmond", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -16,
                "sample", cex = 0.5, facing = "bending")
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -19.5,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -18,
                col = "blanchedalmond", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -18.5,
                "R2", cex = 0.5, facing = "bending")
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -21.8,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -20.3,
                col = "blanchedalmond", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -21,
                "SCD", cex = 0.5, facing = "bending")
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -24.6,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -22.3,
                col = "blanchedalmond", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -23.6,
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


col_response <- structure(c("#15CE59", "#95641A", "#CB3BBD", "#E6880D", "#000000", "#EE0F16", "#8E909B"), names = unique(RF_response$response))
circos.heatmap(RF_response, col = col_response, track.height = 0.05)


### legends
library(ComplexHeatmap)
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

h = dev.size()[2]
lgd_list = packLegend(lgd_Imp, lgd_dir, lgd_R, lgd_SCD,
                      lgd_resp,lgd_sample, max_height = unit(0.9*h, "inch"))
draw(lgd_list, x = circle_size, just = "left")


#Create the circos links
RF_Imp_direction <- RF_Imp.mt
RF_Imp_direction[,c(2,4)] <- NULL
RF_Imp_direction$predictor <- ""
RF_Imp_direction$direction <- ""
RF_Imp_direction$top_predictor <- ""
for (id in unique(RF_Imp_direction$ID)){
  top_pred <- names(which.max(RF_Imp_direction[RF_Imp_direction$ID == id,5:10]))
  RF_Imp_direction[RF_Imp_direction$ID == id,]$top_predictor <- top_pred
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

RF_top_pred <- RF_Imp_direction[,c("ID", "type", "direction", "top_predictor")]
RF_top_pred_name <- RF_top_pred[,c("ID", "top_predictor")]
RF_top_pred_name <- data.frame(RF_top_pred_name[match(RF_Imp.mt.hmap$ID, RF_top_pred_name$ID), ], row.names = 1)

RF_top_pred_dir <- RF_top_pred[,c("ID", "direction")]
RF_top_pred_dir <- data.frame(RF_top_pred_dir[match(RF_Imp.mt.hmap$ID, RF_top_pred_dir$ID), ], row.names = 1)

write.table(RF_Imp_direction, "")

RF_Imp_direction <- filter(RF_Imp_direction, direction != "nonlinear")
RF_Imp_direction$link <- ""

for (s in unique(RF_Imp_direction$sample)) {
  for (p in paste0("H", c(2,4,8,10,11,12))){
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


########### links
# H4 only positive
id <- RF_Imp_direction$ID[RF_Imp_direction$predictor == "H4" & RF_Imp_direction$direction == "pos"]
idx <- which(rownames(RF_Imp) %in% id)
tmp <- t(combn(idx, 2))
df_link_H4 <- data.frame(tmp)
colnames(df_link_H4) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H4))) {
  circos.heatmap.link(df_link_H4$from_index[i],
                      df_link_H4$to_index[i],
                      col = c("#4969D4"), lwd = 0.5)
}

# H2 only positive
id <- RF_Imp_direction$ID[RF_Imp_direction$predictor == "H2" & RF_Imp_direction$direction == "pos"]
idx <- which(rownames(RF_Imp) %in% id)
tmp <- t(combn(idx, 2))
df_link_H2 <- data.frame(tmp)
colnames(df_link_H2) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H2))) {
  circos.heatmap.link(df_link_H2$from_index[i],
                      df_link_H2$to_index[i],
                      col = c("#701717"), lwd = 0.5)
}


# H11 only positive
id <- RF_Imp_direction$ID[RF_Imp_direction$predictor == "H11" & RF_Imp_direction$direction == "pos"]
idx <- which(rownames(RF_Imp) %in% id)
tmp <- t(combn(idx, 2))
df_link_H11 <- data.frame(tmp)
colnames(df_link_H11) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H11))) {
  circos.heatmap.link(df_link_H11$from_index[i],
                      df_link_H11$to_index[i],
                      col = c("#05F3EB"), lwd = 0.5)
}

# H10 only positive
id <- RF_Imp_direction$ID[RF_Imp_direction$predictor == "H10" & RF_Imp_direction$direction == "pos"]
idx <- which(rownames(RF_Imp) %in% id)
tmp <- t(combn(idx, 2))
df_link_H10 <- data.frame(tmp)
colnames(df_link_H10) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H10))) {
  circos.heatmap.link(df_link_H10$from_index[i],
                      df_link_H10$to_index[i],
                      col = c("#71189E"), lwd = 0.5)
}

# H12 only positive
id <- RF_Imp_direction$ID[RF_Imp_direction$predictor == "H12" & RF_Imp_direction$direction == "pos"]
idx <- which(rownames(RF_Imp) %in% id)
tmp <- t(combn(idx, 2))
df_link_H12 <- data.frame(tmp)
colnames(df_link_H12) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H12))) {
  circos.heatmap.link(df_link_H12$from_index[i],
                      df_link_H12$to_index[i],
                      col = c("#890269"), lwd= 0.5)
}

# H8 only positive
id <- RF_Imp_direction$ID[RF_Imp_direction$predictor == "H8" & RF_Imp_direction$direction == "pos"]
idx <- which(rownames(RF_Imp) %in% id)
tmp <- t(combn(idx, 2))
df_link_H8 <- data.frame(tmp)
colnames(df_link_H8) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H8))) {
  circos.heatmap.link(df_link_H8$from_index[i],
                      df_link_H8$to_index[i],
                      col = c("#132892"), lwd = 0.5)
}

write.table(RF_Imp_direction, "")


################# General circos with boxplots ###################

# get the approporitate importance fractions dataframe
library(reshape2)
RF_Imp.t <- RF_Imp.mt.hmap[,c(1:7, 9, 11, 12)]
RF_Imp.t$response <- factor(RF_Imp.t$response, levels = c("H1", "H13", "H3", "H5", "H6", "H7", "H9"))
RF_Imp.t <- arrange(RF_Imp.t, response, type)
RF_Imp.t <- melt(RF_Imp.t[,c(1:7)], id.vars = "ID")
RF_Imp.t <- merge(RF_Imp.t, RF_Imp.mt.hmap[,c("ID", "type", "sample", "response")], by= "ID")
RF_Imp.t$ID <- paste0(RF_Imp.t$sample, "_", RF_Imp.t$variable)
RF_Imp.t$sample <- NULL

# wide
library(tidyr)
RF_Imp.t.wide <- spread(RF_Imp.t, response, value)

RF_Imp.t.wide <- arrange(RF_Imp.t.wide, variable, type)

RF_Imp.t.wide <- RF_Imp.t.wide[,-1] %>% group_by(type, variable) %>% summarise_all(mean)


# get the approporitate importance directions dataframe
library(reshape2)
RF_Dir.t <- RF_Imp.mt.hmap[,c(1, 13:18, 9, 11, 12)]
RF_Dir.t$response <- factor(RF_Dir.t$response, levels = c("H1", "H13", "H3", "H5", "H6", "H7", "H9"))
RF_Dir.t <- arrange(RF_Dir.t, response, type)
RF_Dir.t <- melt(RF_Dir.t[,c(1:7)], id.vars = "ID")
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

RF_Imp.t <- data.matrix(RF_Imp.t.wide[,c(3:9)])
RF_Dir.t <- data.matrix(RF_Dir.t.wide[,c(3:9)])                        

circos.clear()

circos.par(start.degree = 42, gap.after = c(3,3,3,3,3,3,3,3,30,2))

col_fun2 = colorRamp2(c(-1, 0, 1), c("blue", "white", "darkgreen"))
circos.heatmap(RF_Dir.t, split = split, col = col_fun2, track.height = 0.2,
               bg.border = "darkred", bg.lwd = 2, bg.lty = 3, show.sector.labels = TRUE, cluster=T,
               dend.side = "outside", clustering.method = "average", distance.method = "canberra")

circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 9) { # the last sector
    # predictors left
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -3.2,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -1.3,
                col = "#15CE59", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -5.1,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -3.2,
                col = "#95641A", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -7,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -5.1,
                col = "#CB3BBD", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -8.9,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -7,
                col = "#E6880D", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -10.8,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -8.9,
                col = "#000000", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -12.7,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -10.8,
                col = "#EE0F16", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -14.6,
                CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -12.7,
                col = "#8E909B", border = NA)
    # importance annotation
    circos.rect(CELL_META$cell.xlim[2] + convert_x(5.25, "mm"), -14.6,
                CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -1.3,
                col = "blanchedalmond", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -3,
                "TME", cex = 0.4, facing = "bending")
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -6,
                "Response", cex = 0.4, facing = "bending")
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -9,
                "Dependency", cex = 0.38, facing = "bending")
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -12,
                "Direction", cex = 0.38, facing = "bending")
    
    # predictors right
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -3.2,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -1.3,
                col = "#15CE59", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -5.1,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -3.2,
                col = "#95641A", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -7,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -5.1,
                col = "#CB3BBD", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -8.9,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -7,
                col = "#E6880D", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -10.8,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -8.9,
                col = "#000000", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -12.7,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -10.8,
                col = "#EE0F16", border = NA)
    circos.rect(CELL_META$cell.xlim[2] + convert_x(16.75, "mm"), -14.6,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -12.7,
                col = "#8E909B", border = NA)
    
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -18.2,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -15.8,
                col = "blanchedalmond", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -16.7,
                "predictor", cex = 0.5, facing = "bending")
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), -33,
                CELL_META$cell.xlim[2] + convert_x(21, "mm"), -19,
                col = "blanchedalmond", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -21,
                "Mean", cex = 0.5, facing = "bending")
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -24,
                "Fraction", cex = 0.5, facing = "bending")
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -27,
                "Feature", cex = 0.5, facing = "bending")
    circos.text(CELL_META$cell.xlim[2] + convert_x(11, "mm"), -30,
                "Importance", cex = 0.45, facing = "bending")
  }
}, bg.border = NA)



col_pred <- structure(c("#71189E", "#05F3EB", "#890269", "#701717", "#4969D4", "#132892"), names = unique(as.character(RF_Dir.t.wide$variable)))
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
avg_pred <- RF_Imp.mt.hmap[,c(2:7, 11)]
avg_pred <- avg_pred %>% group_by(type) %>% summarise_all(mean)
avg_pred$top_pred <- ""

for(r in 1:nrow(avg_pred)){
  avg_pred[r,]$top_pred<- names(which.max(avg_pred[r,-1]))
}

avg_pred$variable <- paste0(avg_pred$top_pred, "_cor")

RF_Dir.t.wide <- merge(RF_Dir.t.wide, avg_pred[,c("type", "top_pred", "variable")], by = c("type", "variable"), all = T)
colnames(RF_Dir.t.wide)[10] <- "link"

### Add links for direction
RF_Dir.t.wide$link <- ""

for (t in unique(RF_Dir.t.wide$type)){
  for (p in unique(RF_Dir.t.wide$variable)){
    if (rowSums(filter(RF_Dir.t.wide, type == t & variable == p)[,c(3:9)] >= 0.2) >= 4){
      RF_Dir.t.wide[RF_Dir.t.wide$type == t & RF_Dir.t.wide$variable == p,]$link <- paste0(p, "_link")
    }
  }
}

table(RF_Dir.t.wide$link)

# H2 link
id <- rownames(RF_Dir.t.wide)[RF_Dir.t.wide$link == "H2_cor_link" & !is.na(RF_Dir.t.wide$link)]
idx <- which(rownames(as.data.frame(RF_Dir.t)) %in% id)
tmp <- t(combn(idx, 2))
df_link_H2 <- data.frame(tmp)
colnames(df_link_H2) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H2))) {
  circos.heatmap.link(df_link_H2$from_index[i],
                      df_link_H2$to_index[i],
                      col = c("#701717"), lwd= 3)
}

# H11 link
id <- rownames(RF_Dir.t.wide)[RF_Dir.t.wide$link == "H11_cor_link" & !is.na(RF_Dir.t.wide$link)]
idx <- which(rownames(as.data.frame(RF_Dir.t)) %in% id)
tmp <- t(combn(idx, 2))
df_link_H11 <- data.frame(tmp)
colnames(df_link_H11) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H11))) {
  circos.heatmap.link(df_link_H11$from_index[i],
                      df_link_H11$to_index[i],
                      col = c("#05F3EB"), lwd= 3)
}

# H10 link
id <- rownames(RF_Dir.t.wide)[RF_Dir.t.wide$link == "H10_cor_link" & !is.na(RF_Dir.t.wide$link)]
idx <- which(rownames(as.data.frame(RF_Dir.t)) %in% id)
tmp <- t(combn(idx, 2))
df_link_H10 <- data.frame(tmp)
colnames(df_link_H10) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H10))) {
  circos.heatmap.link(df_link_H10$from_index[i],
                      df_link_H10$to_index[i],
                      col = c("#71189E"), lwd= 3)
}




# H4 link
id <- rownames(RF_Dir.t.wide)[RF_Dir.t.wide$link == "H4_cor_link" & !is.na(RF_Dir.t.wide$link)]
idx <- which(rownames(as.data.frame(RF_Dir.t)) %in% id)
tmp <- t(combn(idx, 2))
df_link_H4 <- data.frame(tmp)
colnames(df_link_H4) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H4))) {
  circos.heatmap.link(df_link_H4$from_index[i],
                      df_link_H4$to_index[i],
                      col = c("#4969D4"), lwd= 3)
}


# H8 link
id <- rownames(RF_Dir.t.wide)[RF_Dir.t.wide$link == "H8_cor_link" & !is.na(RF_Dir.t.wide$link)]
idx <- which(rownames(as.data.frame(RF_Dir.t)) %in% id)
tmp <- t(combn(idx, 2))
df_link_H8 <- data.frame(tmp)
colnames(df_link_H8) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H8))) {
  circos.heatmap.link(df_link_H8$from_index[i],
                      df_link_H8$to_index[i],
                      col = c("#132892"), lwd= 3)
}

# H12 link
id <- rownames(RF_Dir.t.wide)[RF_Dir.t.wide$link == "H12_cor_link" & !is.na(RF_Dir.t.wide$link)]
idx <- which(rownames(as.data.frame(RF_Dir.t)) %in% id)
tmp <- t(combn(idx, 2))
df_link_H12 <- data.frame(tmp)
colnames(df_link_H12) <- c("from_index", "to_index")

for(i in seq_len(nrow(df_link_H12))) {
  circos.heatmap.link(df_link_H12$from_index[i],
                      df_link_H12$to_index[i],
                      col = c("#890269"), lwd= 3)
}