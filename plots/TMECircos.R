##################################################
## Project: Cancer Hallmarks
## Script purpose: Plot results of Random Forest models for TME Hallmarks
## Author: Sergi Cervilla* & Mustafa Sibai*
##################################################

library(circlize)
library(paletteer)
library(dplyr)
library(treeshap)
library(stringr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(gridBase)
library(UpSetR)
library(scales)

source("../utils/PlottingMod.R")
source("../utils/SamplesMetadata.R")

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
RF_results_TME <- sapply(RF_Imp.mt.t, function(x) {
  x = x / sum(x)
})
#set rownames
rownames(RF_results_TME) <- rownames(RF_Imp.mt.t)
#transpose again and convert to data frame
RF_results_TME <- data.frame(t(RF_results_TME))
#set rownames as a new column (ID)
RF_results_TME <- tibble::rownames_to_column(RF_results_TME, "ID")
#merge importance fraction and metadata and shapley dependencies
RF_results_TME <- merge(RF_results_TME, RF_Imp.mt[c(1:6, 13:18)], by= "ID")
RF_results_TME$type <- factor(RF_results_TME$type, levels = c("Colorectal", "Liver", "Pancreas", "Breast", "Ovary",  "Kidney", "Bladder", "Glioblastoma", "Prostate", "Lung"))
RF_results_TME <- arrange(RF_results_TME, response,type)

#save tumor types names that will split the heatmap
split <- RF_results_TME$type
split <- factor(split, levels = unique(RF_results_TME$type))

#create a matrix for each information layer of the circos
RF_Imp <- data.matrix(data.frame(RF_results_TME[,c(1:7)], row.names = 1))
RF_Cor <- data.matrix(data.frame(RF_results_TME[,c(1,13:18)], row.names = 1))
RF_SCD <- data.matrix(data.frame(RF_results_TME[,c(1,10)], row.names = 1))
RF_response <- data.frame(RF_results_TME[,c(1,12)], row.names = 1)
RF_R <- data.frame(RF_results_TME[,c(1,8)], row.names = 1)
RF_sample <- data.frame(RF_results_TME[,c(1,9,11,12)], row.names = 1)
RF_sample <- arrange(RF_sample, type, response, sample)

#number the samples within each tumor type
sample_ID <- read.table("../utils/sample_ID.txt", header = T)
RF_sample <- merge(RF_sample, sample_ID, by="sample")
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

lgd_Imp = Legend(title = "Fraction of Importances", col_fun = col_fun1)
lgd_dir = Legend(title = "Feature Dependency", col_fun =  col_fun2)
lgd_R = Legend(title = "R2", col_fun = col_fun3)
lgd_SCD = Legend(title = "Spatial Continuity Degree", col_fun = col_fun4)
lgd_resp = Legend(title = "Response", at = names(col_response),
                  legend_gp = gpar(fill = col_response))
lgd_sample = Legend(title = "Sample", at = names(col_sample),
                    legend_gp = gpar(fill = col_sample))


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



#dendrogram
TME_input_dend <- RI_hallmarks_TME[, c(2,4,5,6,8)]
TME_input_dend <- spread(TME_input_dend, predictor, RI)
TME_input_dend[1:2] <- NULL
TME_input_dend <- TME_input_dend %>% group_by(type) %>%
  dplyr::summarize_all(funs(mean(., na.rm= TRUE)))

TME_input_dend <- data.frame(TME_input_dend, row.names = 1)


# Distance matrix
d <- dist(TME_input_dend, method = "canberra")

# Hierarchical clustering dendrogram
dend <- as.dendrogram(hclust(d)) 

spot_info <- read.table("")
n_type <- spot_info %>% select(sample) %>% unique()
n_type$tumor <- sapply(n_type$sample, function(x) annotation_tissue[[x]])
n_type <- n_type %>% group_by(tumor) %>% summarise(n=n())



labels = TRUE
labels_track_height = 0.1
dend_track_height = 0.5



n_labels <- nleaves(dend)
dend_labels_colors <- labels_colors(dend)
dend_labels_cex <- labels_cex(dend)
labels_dend <- labels(dend)
dend_labels_cex <- rep(1, n_labels)

xl <- numeric(10)
for (i in 1:9) {
  xl[i+1] <- xl[i] + as.numeric(n_type[n_type$tumor == labels_dend[i], "n"]) 
  
}

for (i in 1:10) {
  xl[i] <- xl[i] + as.numeric(n_type[n_type$tumor == labels_dend[i], "n"])/2 
}

pdf("", width = 8, height = 5)
circos.clear()
circos.par(gap.degree=30, start.degree=40)
circlize::circos.initialize("dendrogram", xlim = c(0, sum(n_type$n)))

if (labels) {
  circlize::circos.track(ylim = c(0, 1), panel.fun = function(x, 
                                                              y) {
    circlize::circos.text(x = xl, y = rep(0, 
                                          n_labaels), labels = labels_dend, cex = dend_labels_cex, 
                          col = dend_labels_colors, facing = "clockwise", 
                          niceFacing = TRUE, adj = c(0, 0.5))
  }, bg.border = NA, track.height = labels_track_height)
}

max_height <- attr(dend, "height")
circlize::circos.track(ylim = c(0, max_height), panel.fun = function(x, 
                                                                     y) {
  circos.dendrogram2(dend, facing = facing, max_height = max_height, xl=xl)
}, track.height = dend_track_height, bg.border = NA)

circlize::circos.clear()
invisible(dend)
dev.off()



