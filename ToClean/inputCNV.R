##################################################
## Project: Cancer Hallmarks
## Script purpose: Script to create CNV experiment
## Date: 22/12/2022
## Author: Sergi Cervilla & Mustafa Sibai
##################################################

#All the samples
files <- list.files("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/output_df", full.names = F)
files <- stringr::str_remove(files, pattern = ".txt")
#Exctract count data from the filtered sample (Spot Level)
for (sample in files) {
  STobject <- readRDS(paste0("Desktop/IJC/datasets/Public/", sample, "/RDS/", sample, ".rds"))
  mat <- as.data.frame(STobject@assays[[1]]@counts)
  colnames(mat) <- paste0(STobject@images[[1]]@coordinates[,"row"], "x",STobject@images[[1]]@coordinates[,"col"])
  dir.create(paste0("Desktop/CNV/out/", sample))
  write.table(mat,paste0("Desktop/CNV/input/", sample, ".txt"), quote = F, sep = "\t")
  gc()
}

##############
#Run STARCH (python script)
############

#Add CNV clusters to enhanced Visium sub-spots
for (sample in files) {
  print(sample)
  cnv_clusters <- read.csv(paste0("Desktop/CNV/out/", sample, "/labels_name.csv"), row.names = 1)
  sce <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/enhanced/", sample, "_sce_enhanced.rds"))
  meta <- as.data.frame(colData(sce))  
  rm(sce)
  gc()
  for (spot.id in rownames(cnv_clusters)) {
    coord <- as.integer(strsplit(spot.id, split = "x")[[1]])
    spots <- rownames(meta[meta$spot.row==coord[1] & meta$spot.col==coord[2],])
    df_all[df_all$sample == sample & df_all$spotid %in% spots, "cnv_cluster"] <- cnv_clusters[spot.id,]
  }
}   
write.table(df_all, "Desktop/df_cnv.txt", sep = "\t")




## CNV analysis

#metadata for samples
source("Desktop/IJC/scripts/HallmarksPaper/utils/SamplesMetadata.R")
df_all <- read.table("Desktop/df_cnv.txt")
cnv <- df_all[!is.na(df_all$cnv_cluster), ]


library(tidyverse)
library(tidyr)
#Classify spots in TME, Cancer and Buffer depending on the estimate cluster
cnv$annot <- "TME"
cnv$annot[cnv$clusters == 3] <- "Buffer"
cnv$annot[cnv$clusters %in% c(1,2)] <- "Cancer"
#compute the percentage of TME, Buffer and Cancer spots in each of the CNV clusters
perc <- cnv %>% group_by(sample, cnv_cluster, annot) %>% 
  count() %>% group_by(sample, cnv_cluster) %>% mutate(percent=n/sum(n))  %>% 
  select(-n) %>% spread(annot, value = percent) %>% 
  mutate(BufferCancer = Buffer+Cancer) %>% filter(!is.na(cnv_cluster))
#select all CNV cluster with high content of Cancer (and Buffer) spots
perc_filter <- perc[(perc$Cancer > 0.8 | (perc$BufferCancer > 0.65) & perc$Cancer > 0.35),] %>% filter(!is.na(sample))
#select all spots and clusters that fulfill the above conditions
sub_cnv <- data.frame()
for (sample_s in unique(perc_filter$sample)) {
  for (cluster in filter(perc_filter, sample == sample_s)) {
    sub_cnv <- rbind(sub_cnv, cnv[cnv$sample==sample_s & cnv$cnv_cluster==cluster & cnv$clusters %in% c(1,2),])
  }
}
sub_cnv <- sub_cnv %>% filter(!is.na(sample))
colnames(sub_cnv) <- colnames(cnv)
sub_cnv$cnv_cluster <- factor(sub_cnv$cnv_cluster)
#compute the maximum activity difference (for each hallmark) across clones
test <- sub_cnv %>% group_by(sample, cnv_cluster) %>% summarise(across(H1:H13, mean)) %>% 
  group_by(sample) %>% summarise(minH1 = min(H1), maxH1=max(H1), diffH1=max(H1)-min(H1),
                                 minH2 = min(H2), maxH2=max(H2), diffH2=max(H2)-min(H2),
                                 minH3 = min(H3), maxH3=max(H3), diffH3=max(H3)-min(H3),
                                 minH4 = min(H4), maxH4=max(H4), diffH4=max(H4)-min(H4),
                                 minH5 = min(H5), maxH5=max(H5), diffH5=max(H5)-min(H5),
                                 minH6 = min(H6), maxH6=max(H6), diffH6=max(H6)-min(H6),
                                 minH7 = min(H7), maxH7=max(H7), diffH7=max(H7)-min(H7),
                                 minH8 = min(H8), maxH8=max(H8), diffH8=max(H8)-min(H8),
                                 minH9 = min(H9), maxH9=max(H9), diffH9=max(H9)-min(H9),
                                 minH10 = min(H10), maxH10=max(H10), diffH10=max(H10)-min(H10),
                                 minH11 = min(H11), maxH11=max(H11), diffH11=max(H11)-min(H11),
                                 minH12 = min(H12), maxH12=max(H12), diffH12=max(H12)-min(H12),
                                 minH13 = min(H13), maxH13=max(H13), diffH13=max(H13)-min(H13))
#create a data frame with these hallmark differences
df_diff <- data.frame(row.names = test$sample, test[, paste0("diffH", c(2,4,8,10,11,12))])

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, max(df_diff)), c("white", "red"))
tissue_annot <- sapply(rownames(df_diff), function(x) {annotation_tissue[[x]]})
col_tissue_map <- setNames(Seurat::DiscretePalette(15, palette = "polychrome") [c(1:3, 5, 6:11)], sort(unique(tissue_annot)))

row_ha = rowAnnotation(`Tumor type` = tissue_annot,
                       col=list(`Tumor type` = col_tissue_map))
palette <- do.call(c,color_codes[c(2,4,8,10,11,12)])
names(palette) <- NULL
top_ha = HeatmapAnnotation(" " = anno_boxplot(df_diff, height = unit(2, "cm"),gp = gpar(fill = palette)))
hm <- Heatmap(df_diff, name = "Difference", border = "black", col=col_fun, show_row_names =T, 
        show_column_names = T,
        left_annotation =  row_ha,
        top_annotation = top_ha, 
        rect_gp = gpar(col = "gray", lwd = 0.1))

draw(hm, heatmap_legend_side = "left")

######
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, max(df_diff)), c("white", "red"))
tissue_annot <- sapply(rownames(df_diff), function(x) {annotation_tissue[[x]]})
col_tissue_map <- setNames(Seurat::DiscretePalette(15, palette = "polychrome") [c(1:3, 5, 6:11)], sort(unique(tissue_annot)))

row_ha = HeatmapAnnotation(`Tumor type` = tissue_annot,
                       col=list(`Tumor type` = col_tissue_map))
palette <- do.call(c,color_codes[c(2,4,8,10,11,12)])
names(palette) <- NULL
top_ha = rowAnnotation(" " = anno_boxplot(df_diff, height = unit(2, "cm"),gp = gpar(fill = palette)))
Heatmap(t(df_diff), name = "Difference", border = "black", col=col_fun, show_row_names =T, 
              show_column_names =T,
              top_annotation =  row_ha,
              left_annotation = top_ha, 
              rect_gp = gpar(col = "gray", lwd = 0.1))



######

library(reshape2)
df_diff.long <- melt(df_diff) 
ggplot(df_diff.long, aes(x=variable, y=value)) + geom_boxplot()

library(ggpubr)
boxplot_H2H11 <- sub_cnv[sub_cnv$sample=="M3" & sub_cnv$cnv_cluster %in% c(0,2), c("H11", "H2", "cnv_cluster")]
boxplot_H2H11 <- 
boxplot_H2H11 <- melt(boxplot_H2H11)
ggpubr::ggboxplot(boxplot_H2H11, x = "cnv_cluster", y = "Hallmark activity", 
                  fill = "") 
table(sub_cnv$cnv_cluster, sub_cnv$clusters) 

ggpubr::ggboxplot(sub_cnv[sub_cnv$sample=="M3",], x = "cnv_cluster", y = "H12", 
          fill = "cnv_cluster") + labs(x="CNV clones", y="Hallmark activity (H12)", fill="CNV clones")
hist_df <- perc_filter %>% group_by(sample) %>% count()
hist(hist_df$n)
ggplot(hist_df, aes(x=n)) + geom_bar() + theme_classic() + labs(x="Number of clones")

#percent 
df_cnv <- data.frame(id = paste0(perc$sample, "_", perc$cnv_cluster), annot = perc$annot, y=perc$percent)
df_cnv <- reshape(df_cnv, idvar="id", timevar="annot", direction="wide")
df_cnv[is.na(df_cnv)] <- 0
rownames(df_cnv) <- df_cnv$id
df_cnv$id <- NULL
annot_df <- data.frame(row.names = rownames(df_cnv), cnv = str_sub(rownames(df_cnv), start = -1))
library(pheatmap)
pheatmap(df_cnv, annotation_row = annot_df, show_rownames = F)

###

#### states experiment

#select samples
threshold <- 0.5
df1 <- test %>% 
  rowwise() %>% 
  mutate(mean_diff=mean(c(diffH2,diffH4,diffH8, diffH10, diffH11,diffH12))) %>%
  filter(mean_diff >= threshold) 
df1$sample #samples to select
df2 <- test %>% 
  rowwise() %>% 
  mutate(mean_diff=mean(c(diffH2,diffH4,diffH8,diffH11,diffH12))) %>%
  filter(mean_diff < threshold & mean_diff > 0) 
df2$sample #samples to select



#select clones
maxmin_clones <- sub_cnv %>% filter(sample %in% df1$sample) %>% group_by(sample, cnv_cluster) %>% summarise(across(H1:H13, mean)) %>% 
  group_by(sample) %>% summarise(minH1 = which.min(H1), maxH1=which.max(H1),
                                 minH2 = which.min(H2), maxH2=which.max(H2),
                                 minH3 = which.min(H3), maxH3=which.max(H3),
                                 minH4 = which.min(H4), maxH4=which.max(H4),
                                 minH5 = which.min(H5), maxH5=which.max(H5),
                                 minH6 = which.min(H6), maxH6=which.max(H6),
                                 minH7 = which.min(H7), maxH7=which.max(H7),
                                 minH8 = which.min(H8), maxH8=which.max(H8),
                                 minH9 = which.min(H9), maxH9=which.max(H9),
                                 minH10 = which.min(H10), maxH10=which.max(H10),
                                 minH11 = which.min(H11), maxH11=which.max(H11),
                                 minH12 = which.min(H12), maxH12=which.max(H12),
                                 minH13 = which.min(H13), maxH13=which.max(H13))
#get most frequent minimum clone
min_clones <- apply(maxmin_clones[,paste0("minH",1:13)], 1, function(x) {
  names(which.max(table(x)))
})
names(min_clones) <- maxmin_clones$sample
#get most frquenet maximum clone
max_clones <- apply(maxmin_clones[,paste0("maxH",1:13)], 1, function(x) {
  names(which.max(table(x)))
})
names(max_clones) <- maxmin_clones$sample


#

states <- read.csv("Desktop/CNV/out/M3/states_name.csv", row.names = 1)
order_genes <- read_tsv("Desktop/CNV/siCNV_GeneOrderFile.tsv", col_names = F)
ref <- read.table("Desktop/CNV/STARCH/hgTables_hg38.txt")
chrom_order <-
  ref %>% arrange(factor(chrom))

ref <- ref[order(ref$chrom, ref$cdsStart),]
ref <- ref[!str_detect(ref$chrom, pattern = "_"),]
ref <- unique(ref[,1:2])


order_genes <- ref$name2[ref$name2 %in% rownames(states)]
states <- states[order_genes,]
chr <- data.frame(row.names = order_genes, chr = ref[ref$name2 %in% rownames(states), "chrom"])
chr
pheatmap(states, scale = "none", cluster_rows = F, annotation_row = chr)

mat_min <- c()
mat_max <- c()
for (sample in maxmin_clones$sample) {
  states <- read.csv(paste0("Desktop/CNV/out/", sample, "/states_name.csv"), row.names = 1)
  clone1 <- perc_filter[perc_filter$sample==sample, "cnv_cluster"][max_clones[sample],] #max clone
  clone2 <- perc_filter[perc_filter$sample==sample, "cnv_cluster"][min_clones[sample],] #max clone

  mat_max <- plyr::rbind.fill(mat_max, data.frame(row.names = sample, t(states[,paste0("X", clone1), drop=F])))
  mat_min <- plyr::rbind.fill(mat_min, data.frame(row.names = sample, t(states[,paste0("X", clone2), drop=F])))
}
rownames(mat_max) <- maxmin_clones$sample
rownames(mat_min) <- maxmin_clones$sample


order_genes <- read_tsv("Desktop/CNV/siCNV_GeneOrderFile.tsv", col_names = F)
ref <- read.table("Desktop/CNV/STARCH/hgTables_hg38.txt")
chrom_order <-
  ref %>% arrange(factor(chrom))

ref <- ref[order(match(ref$chrom, paste0("chr", 1:23)), ref$cdsStart),]
ref <- ref[!str_detect(ref$chrom, pattern = "_"),]
ref <- unique(ref[,1:2])


order_genes <- ref$name2[ref$name2 %in% colnames(mat_max)]
mat_max <- mat_max[,order_genes]
mat_min <- mat_min[,order_genes]
chr <- data.frame(row.names = order_genes, chr = ref[ref$name2 %in% colnames(mat_max), "chrom"])
chr$chr <- factor(chr$chr, levels = paste0("chr", 1:23))
pheatmap(mat_max, scale = "none", cluster_rows = F, cluster_cols = F, annotation_col = chr, show_colnames = F)
pheatmap(mat_min, scale = "none", cluster_rows = F, cluster_cols = F, annotation_col = chr, show_colnames = F)



#### Branch method
rownames(df_diff) <- sapply(rownames(df_diff), function(x){
  sample_names[[x]]
})
pheatmap(df_diff[,paste0("diffH", c(2,4,8,11,12))], scale = "none", color = )
clustering <- hclust(dist(df_diff[,paste0("diffH", c(2,4,8,11,12))]))
plot(clustering)
rect.hclust(tree = clustering, k = 6, )
set.seed(123)
ct <- cutree(clustering, 6)
df_cancer_cnv <- df_diff[,paste0("diffH", c(2,4,8,11,12))]
df_cancer_cnv$cluster <- ct
ggplot(df_cancer_cnv, aes(y=diffH12, x=as.factor(cluster))) + geom_boxplot()

#select samples
highdiff_samples <- names(ct[ct %in% c(1,3,5,6)])

lowdiff_samples <- names(ct[ct == 2])

##### HIGH DIFF
#select clones
maxmin_clones <- sub_cnv %>% filter(sample %in% highdiff_samples) %>% group_by(sample, cnv_cluster) %>% summarise(across(H1:H13, mean)) %>% 
  group_by(sample) %>% summarise(                                 minH2 = which.min(H2), maxH2=which.max(H2),
                                 minH4 = which.min(H4), maxH4=which.max(H4),
                                 minH8 = which.min(H8), maxH8=which.max(H8),
                                 minH11 = which.min(H11), maxH11=which.max(H11),
                                 minH12 = which.min(H12), maxH12=which.max(H12))

#get most frequent minimum clone
min_clones <- apply(maxmin_clones[,paste0("minH",c(2,4,8,11,12))], 1, function(x) {
  names(which.max(table(x)))
})
names(min_clones) <- maxmin_clones$sample
#get most frquenet maximum clone
max_clones <- apply(maxmin_clones[,paste0("maxH",c(2,4,8,11,12))], 1, function(x) {
  names(which.max(table(x)))
})
names(max_clones) <- maxmin_clones$sample


#
mat_min <- c()
mat_max <- c()
for (sample in maxmin_clones$sample) {
  states <- read.csv(paste0("Desktop/CNV/out/", sample, "/states_name.csv"), row.names = 1)
  clone1 <- perc_filter[perc_filter$sample==sample, "cnv_cluster"][max_clones[sample],] #max clone
  clone2 <- perc_filter[perc_filter$sample==sample, "cnv_cluster"][min_clones[sample],] #max clone
  
  mat_max <- plyr::rbind.fill(mat_max, data.frame(row.names = sample, t(states[,paste0("X", clone1), drop=F])))
  mat_min <- plyr::rbind.fill(mat_min, data.frame(row.names = sample, t(states[,paste0("X", clone2), drop=F])))
}
rownames(mat_max) <- maxmin_clones$sample
rownames(mat_min) <- maxmin_clones$sample


order_genes <- read_tsv("Desktop/CNV/siCNV_GeneOrderFile.tsv", col_names = F)
ref <- read.table("Desktop/CNV/STARCH/hgTables_hg38.txt")
chrom_order <-
  ref %>% arrange(factor(chrom))

ref <- ref[order(match(ref$chrom, paste0("chr", 1:23)), ref$cdsStart),]
ref <- ref[!str_detect(ref$chrom, pattern = "_"),]
ref <- unique(ref[,1:2])

order_genes <- ref$name2[ref$name2 %in% colnames(mat_max)]
mat_max <- mat_max[,order_genes]
mat_min <- mat_min[,order_genes]

chr <- data.frame(row.names = order_genes, chr = ref[ref$name2 %in% colnames(mat_max), "chrom"])
chr$chr <- factor(chr$chr, levels = paste0("chr", 1:23))
chr_palette <- DiscretePalette(23)
names(chr_palette) <- paste0("chr", 1:23)
chr_col <- list(chr=chr_palette)

rownames(mat_max) <- sapply(rownames(mat_max), function(x){
  sample_names[[x]]
})
rownames(mat_min) <- sapply(rownames(mat_min), function(x){
  sample_names[[x]]
})

pheatmap(mat_max[order(rownames(mat_max)),], scale = "none", cluster_rows = F, cluster_cols = F, annotation_col = chr, show_colnames = F, annotation_colors = chr_col)
pheatmap(mat_min[order(rownames(mat_min)),], scale = "none", cluster_rows = F, cluster_cols = F, annotation_col = chr, show_colnames = F, annotation_colors = chr_col)


boxplot_df <- data.frame(row.names = rownames(mat_max), duplications_max = apply(mat_max, 1, function(x) {
                          noNA <- x[!is.na(x)]
                          sum(noNA==2)/length(noNA)}),
                          duplications_min =apply(mat_min, 1, function(x) {
                           noNA <- x[!is.na(x)]
                           sum(noNA==2)/length(noNA)}),
                         deletions_max = apply(mat_max, 1, function(x) {
                           noNA <- x[!is.na(x)]
                           sum(noNA==0)/length(noNA)}),
                         deletions_min =apply(mat_min, 1, function(x) {
                           noNA <- x[!is.na(x)]
                           sum(noNA==0)/length(noNA)})
                         )

library(reshape2)
boxplot_df <- melt(boxplot_df[!rownames(boxplot_df) %in% c("Glioblastoma1", "Glioblastoma3", "Glioblastoma5"),])
boxplot_df$clone <- "High"
boxplot_df$CNV <- "Duplication"
boxplot_df$clone[boxplot_df$variable %in% c("duplications_min", "deletions_min")] <- "Low"
boxplot_df$CNV[boxplot_df$variable %in% c("deletions_max", "deletions_min")] <- "Deletion"

library(ggpubr)
library(rstatix)
# Statistical test
stat.test <- boxplot_df %>%
  group_by(CNV) %>%
  t_test(value ~ clone, paired = T) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test
ggboxplot(boxplot_df, x="clone", y="value", fill="clone", facet.by = "CNV") + stat_compare_means(paired = TRUE) +
   labs(x="Clone", fill="Clone", y="% of CNV state")



##### LOW DIFF
#select clones
maxmin_clones <- sub_cnv %>% filter(sample %in% lowdiff_samples) %>% group_by(sample, cnv_cluster) %>% summarise(across(H1:H13, mean)) %>% 
  group_by(sample) %>% summarise(                                 minH2 = which.min(H2), maxH2=which.max(H2),
                                                                  minH4 = which.min(H4), maxH4=which.max(H4),
                                                                  minH8 = which.min(H8), maxH8=which.max(H8),
                                                                  minH11 = which.min(H11), maxH11=which.max(H11),
                                                                  minH12 = which.min(H12), maxH12=which.max(H12))

#get most frequent minimum clone
min_clones <- apply(maxmin_clones[,paste0("minH",c(2,4,8,11,12))], 1, function(x) {
  names(which.max(table(x)))
})
names(min_clones) <- maxmin_clones$sample
#get most frquenet maximum clone
max_clones <- apply(maxmin_clones[,paste0("maxH",c(2,4,8,11,12))], 1, function(x) {
  names(which.max(table(x)))
})
names(max_clones) <- maxmin_clones$sample


#
mat_min <- c()
mat_max <- c()
for (sample in maxmin_clones$sample) {
  states <- read.csv(paste0("Desktop/CNV/out/", sample, "/states_name.csv"), row.names = 1)
  clone1 <- perc_filter[perc_filter$sample==sample, "cnv_cluster"][max_clones[sample],] #max clone
  clone2 <- perc_filter[perc_filter$sample==sample, "cnv_cluster"][min_clones[sample],] #max clone
  
  mat_max <- plyr::rbind.fill(mat_max, data.frame(row.names = sample, t(states[,paste0("X", clone1), drop=F])))
  mat_min <- plyr::rbind.fill(mat_min, data.frame(row.names = sample, t(states[,paste0("X", clone2), drop=F])))
}
rownames(mat_max) <- maxmin_clones$sample
rownames(mat_min) <- maxmin_clones$sample


order_genes <- read_tsv("Desktop/CNV/siCNV_GeneOrderFile.tsv", col_names = F)
ref <- read.table("Desktop/CNV/STARCH/hgTables_hg38.txt")
chrom_order <-
  ref %>% arrange(factor(chrom))

ref <- ref[order(match(ref$chrom, paste0("chr", 1:23)), ref$cdsStart),]
ref <- ref[!str_detect(ref$chrom, pattern = "_"),]
ref <- unique(ref[,1:2])

order_genes <- ref$name2[ref$name2 %in% colnames(mat_max)]
mat_max <- mat_max[,order_genes]
mat_min <- mat_min[,order_genes]

chr <- data.frame(row.names = order_genes, chr = ref[ref$name2 %in% colnames(mat_max), "chrom"])
chr$chr <- factor(chr$chr, levels = paste0("chr", 1:23))
chr_palette <- DiscretePalette(23)
names(chr_palette) <- paste0("chr", 1:23)
chr_col <- list(chr=chr_palette)

rownames(mat_max) <- sapply(rownames(mat_max), function(x){
  sample_names[[x]]
})
rownames(mat_min) <- sapply(rownames(mat_min), function(x){
  sample_names[[x]]
})
pheatmap(mat_max[order(rownames(mat_max)),], scale = "none", cluster_rows = F, cluster_cols = F, annotation_col = chr, show_colnames = F, annotation_colors = chr_col)
pheatmap(mat_min[order(rownames(mat_min)),], scale = "none", cluster_rows = F, cluster_cols = F, annotation_col = chr, show_colnames = F, annotation_colors = chr_col)

boxplot_df <- data.frame(row.names = rownames(mat_max), duplications_max = apply(mat_max, 1, function(x) {
  noNA <- x[!is.na(x)]
  sum(noNA==2)/length(noNA)}),
  duplications_min =apply(mat_min, 1, function(x) {
    noNA <- x[!is.na(x)]
    sum(noNA==2)/length(noNA)}),
  deletions_max = apply(mat_max, 1, function(x) {
    noNA <- x[!is.na(x)]
    sum(noNA==0)/length(noNA)}),
  deletions_min =apply(mat_min, 1, function(x) {
    noNA <- x[!is.na(x)]
    sum(noNA==0)/length(noNA)})
)
library(reshape2)
boxplot_df <- melt(boxplot_df)
boxplot_df$clone <- "High"
boxplot_df$CNV <- "Duplication"
boxplot_df$clone[boxplot_df$variable %in% c("duplications_min", "deletions_min")] <- "Low"
boxplot_df$CNV[boxplot_df$variable %in% c("deletions_max", "deletions_min")] <- "Deletion"

library(ggpubr)
library(rstatix)
# Statistical test
stat.test <- boxplot_df %>%
  group_by(CNV) %>%
  t_test(value ~ clone, paired = T) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test
ggboxplot(boxplot_df, x="clone", y="value", fill="clone", facet.by = "CNV") + stat_compare_means(paired = TRUE) +
  labs(x="Clone", fill="Clone", y="% of CNV state")



#Plot cnv boxplot


samples <- c("BreastA", "M3", "UKF260T", "OV4A", "Ovarian", "HCC1T", "HCC2T", "P3", "Co1", "Co4", "Intestine")
for (sample_id in samples) {
  df <- df_all[df_all$sample == sample_id &  df_all$cnv  %in%  perc_filter$cnv_cluster[perc_filter$sample == sample_id] & df_all$clusters %in% 1:2, ] %>% 
    select(H2,H4,H8,H10,H11,H12, cnv_cluster) %>% melt(id.vars = "cnv_cluster")
  df$variable <- factor(df$variable, levels = c(paste0("H", c(11, 4, 12, 8, 2, 10))))
  
  p <- ggplot(df, aes(x=as.factor(variable), y=value, fill=variable)) + geom_boxplot(outlier.size = 0.25) + 
    geom_hline(aes(yintercept = 0), linetype="dashed", col="gray32") + 
    facet_grid(~as.factor(cnv_cluster)) + 
    theme_classic() + scale_fill_manual(values = c("#05F3EB",  "#4969D4", "#890269", "#132892", "#701717", "#71189E")) + 
    labs(x="", y="Hallmark Activity") + theme(legend.position = "none", axis.text.x = element_blank(),
                                              axis.ticks.x = element_blank(), axis.text.y = element_text(size = 45),
                                              axis.title.y =  element_text(size = 55),
                                              strip.text.x = element_text(size = 0)) 
  # convert to grob
  gp <- ggplotGrob(p) # where p is the original ggplot object
  
  # assign the first 4 right-side facet strips with blue fill
  cnv_colors <- c("#E69F00", "#0072B2", "#009E73")
  gp$heights[6] <- unit(0.5, "cm")
  for (i in 1:(length(perc_filter$cnv_cluster[perc_filter$sample == sample_id]))) {
    grob.i <- grep("strip-t", gp$layout$name)[i]
    gp$grobs[[grob.i]]$grobs[[1]]$children[[1]]$gp$fill <- cnv_colors[i]
  }
  png(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/CNV_fig/boxplot_", sample_id, ".png"), units = "cm", height = 20, width = 30 , res = 220)
  grid::grid.draw(gp)
  dev.off()
  
}

gtable::gtable_show_layout(gp)

##########3
cancer <- read.table("Downloads/RF_direction_cancer.txt")
tme <- read.table("Downloads/RF_Imp_circos_direction.txt")

cancer_samples <- cancer  %>% filter(!is.na(link) & link != "" & response == "H10")  %>% select(sample) %>% unique()
tme_samples <- tme  %>% filter(!is.na(link) & link != "")  %>% select(sample) %>% unique()

df_diff$Links <- "None"
df_diff[cancer_samples$sample, "Links"] <- "CancerLink"
df_diff[tme_samples$sample, "Links"] <- "TMELink"
df_diff[intersect(cancer_samples$sample, tme_samples$sample), "Links"] <- "Both"
df_diff.long <- melt(df_diff)
ggplot(df_diff.long, aes(x=Links, y=value)) + geom_boxplot() + facet_wrap(~variable)

links_annot <- rep("None", 57)
names(links_annot) <- rownames(df_diff)
links_annot[cancer_samples$sample] <- "CancerLink"
links_annot[tme_samples$sample] <- "TMELink"
links_annot[intersect(tme_samples$sample, cancer_samples$sample)] <- "CancerLink"


