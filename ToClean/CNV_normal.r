


sample <- "M3"
ST_sce <- readRDS(paste0("Desktop/IJC/datasets/Public/", sample, "/RDS/", sample, "_sce.rds"))
ST_sce.enhanced <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/", sample, "_sce_enhanced.rds"))

subspot_barcode <- transform(merge(
  x = data.frame(ST_sce@colData)[c("spot", "spot.row", "spot.col")],
  y = cbind(rownames = rownames(data.frame(ST_sce.enhanced@colData)), data.frame(ST_sce.enhanced@colData)),
  by = c("spot.row", "spot.col")), row.names = rownames, rownames = NULL)
subspot_barcode$starch_id <- paste0(subspot_barcode$spot.row, "x", subspot_barcode$spot.col)
subspot_barcode[df_all$spotid[df_all$sample==sample], "estimate_cluster"] <- df_all$clusters[df_all$sample==sample]

subspot_barcode[subspot_barcode$estimate_cluster %in% c(4,5),]

subspot_barcode <- subspot_barcode %>% group_by(starch_id) %>% summarize (estimate_cluster = names(which.max(table(estimate_cluster)))) %>% as.data.frame()
rownames(subspot_barcode) <- paste0("X",subspot_barcode$starch_id)
subspot_barcode <- subspot_barcode[colnames(test),]
which(subspot_barcode$estimate_cluster %in% c(4,5)) - 1

write.csv(which(subspot_barcode$estimate_cluster %in% c(4,5)) - 1, "Desktop/CNV/normal.txt", quote = F, row.names = F, col.names = F)

rown
test <- read.table("Desktop/CNV/input/M3.txt", sep = "\t", row.names = F, quote = F, col.names = F)
test


print(sample)
cnv_clusters <- read.csv(paste0("Desktop/CNV/labels_M3.csv"), row.names = 1)
sce <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/enhanced/", sample, "_sce_enhanced.rds"))
meta <- as.data.frame(colData(sce))  
rm(sce)
gc()
for (spot.id in rownames(cnv_clusters)) {
  coord <- as.integer(strsplit(spot.id, split = "x")[[1]])
  spots <- rownames(meta[meta$spot.row==coord[1] & meta$spot.col==coord[2],])
  df_all[df_all$sample == sample & df_all$spotid %in% spots, "cnv_cluster"] <- cnv_clusters[spot.id,]
}


##
sample <- "M3"
sce <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/enhanced/",sample,"_sce_enhanced.rds"))
sce <- sce[,-1]
sce$clusters <- df_all[df_all$sample==sample, "clusters"]
sce$cnv <- df_all[df_all$sample == sample, "cnv_cluster"]
sce$cnv2 <- df[df$sample == sample, "cnv_cluster"]
v <- .make_triangle_subspots(colData(sce), fill = "cnv")
p1 <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~as.factor(fill))) +  theme_void() + coord_flip() + scale_x_reverse() + scale_y_reverse()
v <- .make_triangle_subspots(colData(sce), fill = "cnv2")
p2 <- ggplot()  + geom_polygon(data=v,  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~as.factor(fill))) +  theme_void() + coord_flip() + scale_x_reverse() + scale_y_reverse()

ref_v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/M3.rds")
hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/M3.rds")
for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}
hires + geom_polygon(data=v,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~as.factor(fill))) +  theme_void() +
  labs(fill="CNV clones")




states <- t(read.csv("Desktop/CNV/out/M3/states_name.csv", row.names = 1))
order_genes <- read_tsv("Desktop/CNV/siCNV_GeneOrderFile.tsv", col_names = F)
ref <- read.table("Desktop/CNV/STARCH/hgTables_hg38.txt")
chrom_order <-
  ref %>% arrange(factor(chrom))

ref <- ref[order(match(ref$chrom, paste0("chr", 1:23)), ref$cdsStart),]
ref <- ref[!str_detect(ref$chrom, pattern = "_"),]
ref <- unique(ref[,1:2])

order_genes <- ref$name2[ref$name2 %in% colnames(states)]
states <- states[,order_genes]

chr <- data.frame(row.names = order_genes, chr = ref[ref$name2 %in% colnames(states), "chrom"])
chr$chr <- factor(chr$chr, levels = paste0("chr", 1:23))
chr_palette <- DiscretePalette(23)
names(chr_palette) <- paste0("chr", 1:23)
chr_col <- list(chr=chr_palette)

pheatmap(states, scale = "none", cluster_cols  = F, annotation_col = chr, show_colnames = F)

states <- t(read.csv("Desktop/CNV/states_M3.csv", row.names = 1))
order_genes <- read_tsv("Desktop/CNV/siCNV_GeneOrderFile.tsv", col_names = F)
ref <- read.table("Desktop/CNV/STARCH/hgTables_hg38.txt")
chrom_order <-
  ref %>% arrange(factor(chrom))

ref <- ref[order(match(ref$chrom, paste0("chr", 1:23)), ref$cdsStart),]
ref <- ref[!str_detect(ref$chrom, pattern = "_"),]
ref <- unique(ref[,1:2])

order_genes <- ref$name2[ref$name2 %in% colnames(states)]
states <- states[,order_genes]

chr <- data.frame(row.names = order_genes, chr = ref[ref$name2 %in% colnames(states), "chrom"])
chr$chr <- factor(chr$chr, levels = paste0("chr", 1:23))
chr_palette <- DiscretePalette(23)
names(chr_palette) <- paste0("chr", 1:23)
chr_col <- list(chr=chr_palette)

pheatmap(states, scale = "none", cluster_cols  = F, annotation_col = chr, show_colnames = F)

states <- t(read.csv("Desktop/CNV/states_M3_1.csv", row.names = 1))
order_genes <- read_tsv("Desktop/CNV/siCNV_GeneOrderFile.tsv", col_names = F)
ref <- read.table("Desktop/CNV/STARCH/hgTables_hg38.txt")
chrom_order <-
  ref %>% arrange(factor(chrom))

ref <- ref[order(match(ref$chrom, paste0("chr", 1:23)), ref$cdsStart),]
ref <- ref[!str_detect(ref$chrom, pattern = "_"),]
ref <- unique(ref[,1:2])

order_genes <- ref$name2[ref$name2 %in% colnames(states)]
states <- states[,order_genes]

chr <- data.frame(row.names = order_genes, chr = ref[ref$name2 %in% colnames(states), "chrom"])
chr$chr <- factor(chr$chr, levels = paste0("chr", 1:23))
chr_palette <- DiscretePalette(23)
names(chr_palette) <- paste0("chr", 1:23)
chr_col <- list(chr=chr_palette)

pheatmap(states, scale = "none", cluster_cols  = F, annotation_col = chr, show_colnames = F)


states <- t(read.csv("Desktop/CNV/states_M3_2.csv", row.names = 1))
order_genes <- read_tsv("Desktop/CNV/siCNV_GeneOrderFile.tsv", col_names = F)
ref <- read.table("Desktop/CNV/STARCH/hgTables_hg38.txt")
chrom_order <-
  ref %>% arrange(factor(chrom))

ref <- ref[order(match(ref$chrom, paste0("chr", 1:23)), ref$cdsStart),]
ref <- ref[!str_detect(ref$chrom, pattern = "_"),]
ref <- unique(ref[,1:2])

order_genes <- ref$name2[ref$name2 %in% colnames(states)]
states <- states[,order_genes]

chr <- data.frame(row.names = order_genes, chr = ref[ref$name2 %in% colnames(states), "chrom"])
chr$chr <- factor(chr$chr, levels = paste0("chr", 1:23))
chr_palette <- DiscretePalette(23)
names(chr_palette) <- paste0("chr", 1:23)
chr_col <- list(chr=chr_palette)

pheatmap(states, scale = "none", cluster_cols  = F, annotation_col = chr, show_colnames = F)
