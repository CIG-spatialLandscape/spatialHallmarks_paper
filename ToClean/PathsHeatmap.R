library(ComplexHeatmap)
library(circlize)

#M3

paths <- read.table("Downloads/M3.txt", sep = "\t")

paths[, 15:ncol(paths)] <- scale(paths[, 15:ncol(paths)])
cnv <- read.table("Desktop/df_cnv.txt")
cnv <- cnv %>% filter(sample == "M3")

paths[cnv$spotid, "cnv_cluster"] <- cnv$cnv_cluster

paths.cnv <- paths %>% filter(estimate.cluster %in% 1:2 & cnv_cluster %in% 0:2) %>%
  select(-estimate.cluster) %>% group_by(cnv_cluster) %>% summarise_all(mean)

mat.paths <- as.data.frame(matrix(nrow = 3, ncol = 0))
hallmark_annot <- c()
for (hallmark in paste0("H", c(2,4,8,10,11,12))) {
  col_paths <- grep(hallmark, colnames(paths.cnv),value = T)
  sub_mat <- paths.cnv[, col_paths[2:length(col_paths)]]
  dist_mat <- dist(t(sub_mat))
  clustering <- hclust(dist_mat)
  mat.paths <- cbind(mat.paths, paths.cnv[, c(hallmark, colnames(sub_mat)[clustering$order])])
  hallmark_annot <- c(hallmark_annot, rep(hallmark, length(col_paths)))
  
}
hallmark_annot <- factor(hallmark_annot, levels = c("H11", "H4", "H12", "H8", "H2", "H10"))
rownames(mat.paths) <- 0:2



palette <- c("#701717",  "#4969D4", "#132892", "#71189E", "#05F3EB", "#890269")
names(palette) <- c("H2", "H4", "H8", "H10", "H11", "H12")
col_ha <- HeatmapAnnotation(Hallmark = hallmark_annot, 
                            col=list(Hallmark=palette))
ha = HeatmapAnnotation(foo = anno_mark(at = which(colnames(mat.paths) %in%  c("H2", "H4", "H8", "H10", "H11", "H12")), 
                                   labels = c("H2", "H4", "H8", "H10", "H11", "H12")))
pdf("Desktop/heatmap_M3.pdf", width = 20, height = 5)
Heatmap(scale(mat.paths), cluster_columns = F, show_column_names = F,
        bottom_annotation = col_ha, top_annotation = ha, column_split = hallmark_annot,
        border = "gray22", column_title_gp = gpar(fontsize = 0
        ), column_gap = unit(0.25, "cm"))

dev.off()

mat_order <- colnames(mat.paths)
#BreastA

paths <- read.table("Downloads/BreastA.txt", sep = "\t")

paths[, 15:ncol(paths)] <- scale(paths[, 15:ncol(paths)])
cnv <- read.table("Desktop/df_cnv.txt")
cnv <- cnv %>% filter(sample == "BreastA")

paths[cnv$spotid, "cnv_cluster"] <- cnv$cnv_cluster

paths.cnv <- paths %>% filter(estimate.cluster %in% 1:2 & cnv_cluster %in% 0:2) %>%
  select(-estimate.cluster) %>% group_by(cnv_cluster) %>% summarise_all(mean)

mat.paths <- as.data.frame(matrix(nrow = 3, ncol = 0))
hallmark_annot <- c()
for (hallmark in paste0("H", c(2,4,8,10,11,12))) {
  col_paths <- grep(hallmark, colnames(paths.cnv),value = T)
  sub_mat <- paths.cnv[, col_paths[2:length(col_paths)]]
  dist_mat <- dist(t(sub_mat))
  clustering <- hclust(dist_mat)
  mat.paths <- cbind(mat.paths, paths.cnv[, c(hallmark, colnames(sub_mat)[clustering$order])])
  hallmark_annot <- c(hallmark_annot, rep(hallmark, length(col_paths)))
}
hallmark_annot <- factor(hallmark_annot, levels = c("H11", "H4", "H12", "H8", "H2", "H10"))
rownames(mat.paths) <- 0:2

library(ComplexHeatmap)
library(circlize)


palette <- c("#701717",  "#4969D4", "#132892", "#71189E", "#05F3EB", "#890269")
names(palette) <- c("H2", "H4", "H8", "H10", "H11", "H12")
col_ha <- HeatmapAnnotation(Hallmark = hallmark_annot, 
                            col=list(Hallmark=palette))
ha = HeatmapAnnotation(foo = anno_mark(at = which(colnames(mat.paths) %in%  c("H2", "H4", "H8", "H10", "H11", "H12")), 
                                       labels = c("H2", "H4", "H8", "H10", "H11", "H12")))
pdf("Desktop/heatmap_BreastA.pdf", width = 20, height = 5)
mat.paths <- mat.paths[, mat_order]
Heatmap(scale(mat.paths), cluster_columns = F, show_column_names = F,
        bottom_annotation = col_ha, top_annotation = ha, column_split = hallmark_annot,
        border = "gray22", column_title_gp = gpar(fontsize = 0
        ), column_gap = unit(0.25, "cm"))

dev.off()

#UKF260T

paths <- read.table("Downloads/UKF260T.txt", sep = "\t")

paths[, 15:ncol(paths)] <- scale(paths[, 15:ncol(paths)])
cnv <- read.table("Desktop/df_cnv.txt")
cnv <- cnv %>% filter(sample == "UKF260T")

paths[cnv$spotid, "cnv_cluster"] <- cnv$cnv_cluster

paths.cnv <- paths %>% filter(estimate.cluster %in% 1:2 & cnv_cluster %in% 0:2) %>%
  select(-estimate.cluster) %>% group_by(cnv_cluster) %>% summarise_all(mean)

mat.paths <- as.data.frame(matrix(nrow = 3, ncol = 0))
hallmark_annot <- c()
for (hallmark in paste0("H", c(2,4,8,10,11,12))) {
  col_paths <- grep(hallmark, colnames(paths.cnv),value = T)
  sub_mat <- paths.cnv[, col_paths[2:length(col_paths)]]
  dist_mat <- dist(t(sub_mat))
  clustering <- hclust(dist_mat)
  mat.paths <- cbind(mat.paths, paths.cnv[, c(hallmark, colnames(sub_mat)[clustering$order])])
  hallmark_annot <- c(hallmark_annot, rep(hallmark, length(col_paths)))
}
hallmark_annot <- factor(hallmark_annot, levels = c("H11", "H4", "H12", "H8", "H2", "H10"))
rownames(mat.paths) <- 0:2

library(ComplexHeatmap)
library(circlize)


palette <- c("#701717",  "#4969D4", "#132892", "#71189E", "#05F3EB", "#890269")
names(palette) <- c("H2", "H4", "H8", "H10", "H11", "H12")
col_ha <- HeatmapAnnotation(Hallmark = hallmark_annot, 
                            col=list(Hallmark=palette))
ha = HeatmapAnnotation(foo = anno_mark(at = which(colnames(mat.paths) %in%  c("H2", "H4", "H8", "H10", "H11", "H12")), 
                                       labels = c("H2", "H4", "H8", "H10", "H11", "H12")))
pdf("Desktop/heatmap_UKF260T.pdf", width = 20, height = 5)
mat.paths <- mat.paths[, mat_order]
Heatmap(scale(mat.paths), cluster_columns = F, show_column_names = F,
        bottom_annotation = col_ha, top_annotation = ha, column_split = hallmark_annot,
        border = "gray22", column_title_gp = gpar(fontsize = 0
        ), column_gap = unit(0.25, "cm"))

dev.off()


#Plot M3 pathways

paths <- read.table("Downloads/M3.txt", sep = "\t", check.names = F)

paths[, 15:ncol(paths)] <- scale(paths[, 15:ncol(paths)])
cnv <- read.table("Desktop/df_cnv.txt")
cnv <- cnv %>% filter(sample == "M3")

paths[cnv$spotid, "cnv_cluster"] <- cnv$cnv_cluster

#paths.cnv <- paths %>% filter(estimate.cluster %in% 1:2 & cnv_cluster %in% 0:2) %>%
 # select(-estimate.cluster)
grep("HAT", colnames(paths), value = T)

  pathways <- c("H10_Glycolysis", "H10_Warburg.Effect", "H10_Regulation.of.glycolysis.by.fructose.2.6.bisphosphate.metabolism",
              "H11_Senescence.Associated.Secretory.Phenotype..SASP.", "H11_Oxidative.Stress.Induced.Senescence",
              "H11_Formation.of.Senescence.Associated.Heterochromatin.Foci..SAHF.", "H2_p53.pathway", 
              "H2_SMAD2.SMAD3.SMAD4.heterotrimer.regulates.transcription", "H2_Cdc20.Phospho.APC.C.mediated.degradation.of.Cyclin.A",
              "H8_G2.M.DNA.damage.checkpoint", "H8_Gap.filling.DNA.repair.synthesis.and.ligation.in.GG.NER", "H8_Gap.filling.DNA.repair.synthesis.and.ligation.in.TC.NER",
              "H4_DNA.Damage.Telomere.Stress.Induced.Senescence", "H8_Recognition.of.DNA.damage.by.PCNA.containing.replication.complex",
              "H4_Packaging.Of.Telomere.Ends", "H12_Orc1.removal.from.chromatin","H12_SUMOylation.of.chromatin.organization.proteins",
              "H12_HATs.acetylate.histones"
              )
  
  pathways <- c("H10_Glycolysis", "H11_Senescence-Associated Secretory Phenotype (SASP)",
                "H2_p53 pathway", "H4_Packaging Of Telomere Ends", "H8_G2/M DNA damage checkpoint", "H12_HATs acetylate histones")
  

sample <- "M3"
sce <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/enhanced/", sample, "_sce_enhanced.rds"))
sce <- sce[,-1]
pathway <- pathways[2]
colData(sce)[, c(pathways, "clusters", "cnv_cluster")] <- paths[, c(pathways, "estimate.cluster", "cnv_cluster")]
sce$empty <- "empty"
ref_v <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/",sample,".rds"))
hires <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/",sample,".rds"))
v3 <- .make_triangle_subspots(colData(sce)[!(sce$clusters %in% c(1,2) & sce$cnv_cluster %in% 0:2),], fill = "empty")
for (spot in unique(v3$spot)) {
  v3$imagecol[v3$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v3$imagerow[v3$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}
for (pathway in pathways) {
  print(pathway)
  v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2) & sce$cnv_cluster %in% 0:2,], fill = pathway)
  for (spot in unique(v$spot)) {
    v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
    v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
  }
  hires + geom_polygon(data=v,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
    scale_fill_gradientn("Pathway activity", colours = viridisLite::rocket(1000, alpha = 1, begin = 0, end = 1, direction = 1)[200:1000]) + 
    ggnewscale::new_scale("fill") + geom_polygon(data=v3,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~as.factor(fill))) +
    scale_fill_manual("Buffer", values = "gray41") + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(str_split(pathway, pattern = "_", simplify = T)[1,2])
  ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/CNV_fig/", sample, "_", pathway, ".png"), bg = "white", width = 7, height = 7)
  
}




