library(paletteer)

df1 <- read.table("Downloads/RF_Imp_circos.txt")
df2 <- read.table("Downloads/RF_Imp_circos_direction.txt")
head(df1)


unique(df2$link)
id <- df2$ID[df2$link == "H2_link" & !is.na(df2$link)]
idx <- which(rownames(df1) %in% id)
tmp <- t(combn(idx, 2))
tmp


for (h_cancer in paste0("H", c(2,4,8,11,12))) {
  id <- df2$ID[df2$link == paste0(h_cancer, "_link") & !is.na(df2$link)]
  for (s in id) {
    sample <- str_split(s, pattern = "_")[[1]][1]
    if (file.exists(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/",sample,".rds"))) {
      print(sample)
      h_tme <- str_split(s, pattern = "_")[[1]][2]
      
      sce <- readRDS(paste0("Desktop/enhanced/enhanced_sce/", sample, "_sce_enhanced.rds"))
      sce <- sce[,-1]
      colData(sce)[, c(h_cancer, h_tme, "clusters")] <- df_all[df_all$sample==sample, c(h_cancer, h_tme, "clusters")]
      
      
      v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill = h_cancer)
      v2 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(4,5),], fill = h_tme)
      
      ref_v <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/",sample,".rds"))
      hires <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/",sample,".rds"))
      for (spot in unique(v$spot)) {
        v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
        v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
      }
      
      
      for (spot in unique(v2$spot)) {
        v2$imagecol[v2$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
        v2$imagerow[v2$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
      }
      
      hires + geom_polygon(data=v2,  aes_(x=~imagecol, y=~imagerow, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
        scale_fill_gradientn(h_tme, colours = rev(paletteer_c("grDevices::YlOrRd", 30))  ) + 
        ggnewscale::new_scale("fill") + geom_polygon(data=v,  aes_(x=~imagecol, y=~imagerow, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
        scale_fill_gradientn(h_cancer, colours = rev(
          
          paletteer_c("grDevices::Blues 2", 30)
        ))
      
      ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RF_plots/links/", sample, "_", h_cancer, "_", h_tme, ".png"), bg = "white", width = 7, height = 7)
      
    }
    
  }
}

df1 <- read.table("Downloads/RF_Imp_circos.txt")
df2 <- read.table("Downloads/RF_direction_cancer.txt")
head(df2)



df_all <- read.table("Desktop/df_cnv.txt")
for (h_tme in paste0("H", c(1,3,5,6,7,9,13))) {
  id <- df2$ID[df2$link == paste0(h_tme, "_link") & !is.na(df2$link)]
  for (s in id) {
    sample <- str_split(s, pattern = "_")[[1]][1]
    if (file.exists(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/",sample,".rds"))) {
      print(sample)
      h_cancer <- str_split(s, pattern = "_")[[1]][2]
      
      sce <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/enhanced/", sample, "_sce_enhanced.rds"))
      sce <- sce[,-1]
      colData(sce)[, c(h_cancer, h_tme, "clusters")] <- df_all[df_all$sample==sample, c(h_cancer, h_tme, "clusters")]
      
      
      v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill = h_cancer)
      v2 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(4,5),], fill = h_tme)
      
      ref_v <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/",sample,".rds"))
      hires <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/",sample,".rds"))
      for (spot in unique(v$spot)) {
        v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
        v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
      }
      
      
      for (spot in unique(v2$spot)) {
        v2$imagecol[v2$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
        v2$imagerow[v2$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
      }
      
      hires + geom_polygon(data=v2,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
        scale_fill_gradientn(h_tme, colours = rev(paletteer_c("grDevices::YlOrRd", 30))  ) + 
        ggnewscale::new_scale("fill") + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
        scale_fill_gradientn(h_cancer, colours = rev(
          
          paletteer_c("grDevices::Blues 2", 30)
        ))
      
      ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RF_plots/links/", sample, "_", h_cancer, "_", h_tme, ".png"), bg = "white", width = 7, height = 7)
      
    }
    
  }
}



library(stringr)
for (h_tme in paste0("H", c(1,3,5,6,7,9,13))) {
  id <- df2$ID[df2$link == paste0(h_tme, "_link") & !is.na(df2$link)]
  for (s in id) {
    sample <- str_split(s, pattern = "_")[[1]][1]
    if (!file.exists(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/",sample,".rds")))  print(sample)
  }
}




df_all <- read.table("Desktop/df_cnv.txt")
df_tme <- read.table("Desktop/df_tme_fix.txt")
for (h_tme in paste0("H", c(1,3,5,6,7,9,13))) {
  id <- df2$ID[df2$link == paste0(h_tme, "_link") & !is.na(df2$link)]
  for (s in id) {
    sample <- str_split(s, pattern = "_")[[1]][1]
    if (file.exists(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/",sample,".rds"))) {
      print(sample)
      h_cancer <- str_split(s, pattern = "_")[[1]][2]
      
      sce <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/enhanced/", sample, "_sce_enhanced.rds"))
      sce <- sce[,-1]
      colData(sce)[, c(h_cancer, h_tme, "clusters")] <- df_all[df_all$sample==sample, c(h_cancer, h_tme, "clusters")]
      sce <- sce[,sce$clusters %in% c(4,5)]
      colData(sce)[, paste0(h_cancer, "_cancer")] <- df_tme[df_tme$sample==sample, paste0(h_cancer, "_cancer")]
      
      ref_v <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/",sample,".rds"))
      hires <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/",sample,".rds"))
      
      v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(4,5),], fill = h_tme)
      
      for (spot in unique(v$spot)) {
        v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
        v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
      }
      
      hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
        scale_fill_viridis_c(h_tme,option = "B")
      ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RF_plots/links/tme/", sample, "_", h_tme, ".png"), bg = "white", width = 7, height = 7)
      
      v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(4,5),], fill =  paste0(h_cancer, "_cancer"))
      
      for (spot in unique(v$spot)) {
        v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
        v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
      }
      
      hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
        scale_fill_viridis_c(paste0(h_cancer, "_radar"))
      ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RF_plots/links/tme/", sample, "_", h_cancer, "_radar.png"), bg = "white", width = 7, height = 7)
      
    }
    
  }
}




df_all <- read.table("Desktop/df_cnv.txt")
df_cancer <- read.table("Desktop/df_tme_fix.txt")
for (h_tme in paste0("H",  c(2,4,8,11,12))) {
  id <- df2$ID[df2$link == paste0(h_cancer, "_link") & !is.na(df2$link)]
  for (s in id) {
    sample <- str_split(s, pattern = "_")[[1]][1]
    if (file.exists(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/",sample,".rds"))) {
      print(sample)
      h_tme <- str_split(s, pattern = "_")[[1]][2]
      
      sce <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/enhanced/", sample, "_sce_enhanced.rds"))
      sce <- sce[,-1]
      colData(sce)[, c(h_cancer, h_cancer, "clusters")] <- df_all[df_all$sample==sample, c(h_cancer, h_tme, "clusters")]
      sce <- sce[,sce$clusters %in% c(1,2)]
      colData(sce)[, paste0(h_TME, "_TME")] <- df_tme[df_tme$sample==sample, paste0(h_tme, "_TME")]
      
      ref_v <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/",sample,".rds"))
      hires <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/",sample,".rds"))
      
      v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill = h_tme)
      
      for (spot in unique(v$spot)) {
        v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
        v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
      }
      
      hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
        scale_fill_viridis_c(h_cancer,option = "B")
      ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RF_plots/links/tme/", sample, "_", h_tme, ".png"), bg = "white", width = 7, height = 7)
      
      v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill =  paste0(h_cancer, "_cancer"))
      
      for (spot in unique(v$spot)) {
        v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
        v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
      }
      
      hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
        scale_fill_viridis_c(paste0(h_tme, "_radar"))
      ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RF_plots/links/tme/", sample, "_", h_cancer, "_radar.png"), bg = "white", width = 7, height = 7)
      
    }
  }
}




################################################################################

df_all <- read.table("Desktop/df_cnv.txt")
for (h_tme in paste0("H", c(1,3,5,6,7,9,13))) {
  id <- df2$ID[df2$link == paste0(h_tme, "_link") & !is.na(df2$link)]
  for (s in id) {
    sample <- str_split(s, pattern = "_")[[1]][1]
    if (file.exists(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/",sample,".rds"))) {
      
      
      
      
      print(sample)
      h_cancer <- str_split(s, pattern = "_")[[1]][2]
      
      sce <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/enhanced/", sample, "_sce_enhanced.rds"))
      sce <- sce[,-1]
      colData(sce)[, c(h_cancer, h_tme, "clusters")] <- df_all[df_all$sample==sample, c(h_cancer, h_tme, "clusters")]
      colData(sce)[, c("H12")] <- df_all[df_all$sample==sample, c("H12")]
      sce$diff <- sce$H10 - sce$H8
      sce$diff <- sce$H12 - sce$H8
      
      v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill = "diff")
      v2 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(4,5),], fill = h_tme)
      
      ref_v <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/",sample,".rds"))
      hires <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/",sample,".rds"))
      for (spot in unique(v$spot)) {
        v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
        v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
      }
      
      
      for (spot in unique(v2$spot)) {
        v2$imagecol[v2$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
        v2$imagerow[v2$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
      }
      lbl <- c("H12"=max(v2$fill)-0.5, "0"=0, "H8"=min(v2$fill)+0.5)
      
      
      hires + geom_polygon(data=v2,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
        scale_fill_viridis_c(h_tme, option = "B") + 
        ggnewscale::new_scale("fill") + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
        scale_fill_gradient2("diff",low = "darkmagenta", mid = "white", high = "forestgreen", breaks  = lbl)

      
      
      
      
      
      hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
        scale_fill_gradientn("diff", colours = rev(
          
          paletteer_c("grDevices::Blues 2", 30)
        ))
      
      ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RF_plots/links/", sample, "_", h_cancer, "_", h_tme, ".png"), bg = "white", width = 7, height = 7)
      
    }
    
  }
}


sample <- "P306"
h_cancer1 <- "H8"
h_cancer2 <- "H10"
h_tme <- "H1"

sce <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/enhanced/", sample, "_sce_enhanced.rds"))
sce <- sce[,-1]
colData(sce)[, c(h_cancer1, h_cancer2, h_tme, "clusters")] <- df_all[df_all$sample==sample, c(h_cancer1, h_cancer2, h_tme, "clusters")]
sce$diff <- colData(sce)[, h_cancer1] - colData(sce)[, h_cancer2] 


v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill = "diff")
v2 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(4,5),], fill = h_tme)
v3 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(3),], fill = "clusters")

ref_v <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/",sample,".rds"))
hires <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/",sample,".rds"))

for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}


for (spot in unique(v2$spot)) {
  v2$imagecol[v2$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v2$imagerow[v2$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}

for (spot in unique(v3$spot)) {
  v3$imagecol[v3$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v3$imagerow[v3$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}

hires + geom_polygon(data=v2,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn(h_tme, colours = viridisLite::rocket(1000, alpha = 1, begin = 0, end = 1, direction = 1)[200:1000]) + 
  ggnewscale::new_scale("fill") + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +
  scale_fill_gradient2(paste0(h_cancer1, "-", h_cancer2), high = "green4", mid = "white", low  = "navyblue", limits = c(min(sce$diff),max(sce$diff))) + 
   ggnewscale::new_scale("fill") + geom_polygon(data=v3,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~as.factor(fill))) +
  scale_fill_manual("Buffer", values = "gray41") 

ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RF_plots/examples/", sample, "_", h_tme, "_", h_cancer1, h_cancer2,".png"), bg = "white", width = 7, height = 7)

#####

rowSums(table(df2$sample, df2$link)[,-1])



df2[df2$link != "" & !is.na(df2$link),]


df_hallmarks <- df_all[, c(paste0("H", 1:13),"clusters", "sample")]
df_hallmarks$compartment <- cut(df_hallmarks$clusters, breaks = c(0,2.5,3.5,6), labels = c("Cancer", "Buffer", "TME"))

df2 <- read.table("Downloads/RF_direction_cancer.txt")
#df2 <- read.table("Downloads/RF_Imp_circos_direction.txt")
df_hallmarks.long <- df_hallmarks %>% select(-clusters) %>% group_by(sample, compartment) %>% 
  summarise_all(mean) %>% ungroup() %>% melt(id.vars = c("sample", "compartment")) %>%
  mutate(ID = paste0(sample, "_", variable)) %>% select(ID, compartment, value)
df_links <- merge(df2, df_hallmarks.long[df_hallmarks.long$compartment == "Cancer", c("ID", "value")], by ="ID")


colnames(df_hallmarks)[1:13] <- hallmark_names
df_hallmarks.long <- melt(df_hallmarks, id.vars = c("sample", "clusters"))
df_hallmarks.long$ID <- paste0(df_hallmarks.long$sample, "_", df_hallmarks.long)


####
df2 <- read.table("Downloads/RF_Imp_direction_cancer.txt")
df2 <- read.table("Downloads/RF_Imp_direction_TME.txt")

sample <- "HCC2T"
h_tme1 <- "H1"
h_tme2 <- "H6"
h_cancer <- "H12"

sce <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/enhanced/", sample, "_sce_enhanced.rds"))
sce <- sce[,-1]
colData(sce)[, c(h_tme1, h_tme2, h_cancer, "clusters")] <- df_all[df_all$sample==sample, c(h_tme1, h_tme2, h_cancer, "clusters")]
sce$diff <- colData(sce)[, h_tme1] - colData(sce)[, h_tme2] 


v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(4,5),], fill = "diff")
v2 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill = h_cancer)
v3 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(3),], fill = "clusters")

ref_v <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/",sample,".rds"))
hires <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/",sample,".rds"))

for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}


for (spot in unique(v2$spot)) {
  v2$imagecol[v2$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v2$imagerow[v2$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}

for (spot in unique(v3$spot)) {
  v3$imagecol[v3$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v3$imagerow[v3$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}

hires + geom_polygon(data=v2,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn(h_cancer, colours = viridisLite::rocket(1000, alpha = 1, begin = 0, end = 1, direction = 1)[200:1000]) + 
  ggnewscale::new_scale("fill") + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +
  scale_fill_gradient2(paste0(h_tme1, "-", h_tme2), high = "green4", mid = "white", low  = "navyblue", limits = c(min(sce$diff),max(sce$diff))) + 
  ggnewscale::new_scale("fill") + geom_polygon(data=v3,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~as.factor(fill))) +
  scale_fill_manual("Buffer", values = "gray41") 

ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RF_plots/examples/", sample, "_", h_cancer, "_", h_tme1, h_tme2,".png"), bg = "white", width = 7, height = 7)
#####

df_tme_direction <- read.table("Downloads/RF_Imp_direction_TME.txt")

#####
samples <- c("DU13", "P4", "OVFFPE", "HCC2T", "M4", "Ovarian", "Glioblastoma", 
             "PC1", "PC2", "P288", "P306", "M1", "BreastA", "C20", "Co1", "Co4")
for (sample in samples) {
  print(sample)
  sce <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/enhanced/", sample, "_sce_enhanced.rds"))
  sce <- sce[,-1]
  colData(sce)[, "clusters"] <- df_all[df_all$sample==sample, "clusters"]
  sce$compartments <- cut(sce$clusters, breaks = c(0,2.5,3.5,5), labels = c("Cancer", "Buffer", "TME"))
  v <- .make_triangle_subspots(colData(sce), fill = "compartments")
  ref_v <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/",sample,".rds"))
  hires <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/",sample,".rds"))
  for (spot in unique(v$spot)) {
    v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
    v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
  }
  hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~factor(fill))) +  theme_void() + coord_equal() +
    scale_fill_manual(values = rev(c("lightgoldenrod1","lightpink2", "orchid4"))) + labs(fill="Compartments")  
  ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RF_plots/estimate_clusters/", sample, ".png"), bg = "white", width = 7, height = 7)
}


###



####
df2 <- read.table("Downloads/RF_Imp_direction_cancer.txt")
df2 <- read.table("Downloads/RF_Imp_direction_TME.txt")

sample <- "HCC2T"
h_tme1 <- "H6"
h_tme2 <- "H6"
h_cancer <- "H2"

sce <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/enhanced/", sample, "_sce_enhanced.rds"))
sce <- sce[,-1]
colData(sce)[, c(h_tme1, h_tme2, h_cancer, "clusters")] <- df_all[df_all$sample==sample, c(h_tme1, h_tme2, h_cancer, "clusters")]
sce$diff <- colData(sce)[, h_tme1] - colData(sce)[, h_tme2] 


v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(4,5),], fill = h_tme1)
v2 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill = h_cancer)
v3 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(3),], fill = "clusters")

ref_v <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/",sample,".rds"))
hires <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/",sample,".rds"))

for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}


for (spot in unique(v2$spot)) {
  v2$imagecol[v2$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v2$imagerow[v2$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}

for (spot in unique(v3$spot)) {
  v3$imagecol[v3$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v3$imagerow[v3$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}

hires + geom_polygon(data=v2,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn(h_cancer, colours = viridisLite::rocket(1000, alpha = 1, begin = 0, end = 1, direction = 1)[200:1000]) + 
  ggnewscale::new_scale("fill") + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +
  scale_fill_gradientn(h_tme1, colours = viridisLite::mako(1000, alpha = 1, begin = 0, end = 1, direction = 1)[200:1000]) + 
  ggnewscale::new_scale("fill") + geom_polygon(data=v3,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~as.factor(fill))) +
  scale_fill_manual("Buffer", values = "gray41") 


ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RF_plots/examples/", sample, "_", h_cancer, "_", h_tme1,".png"), bg = "white", width = 7, height = 7)


####

df_all <- read.table("Desktop/df_cnv.txt")
sample <- "M3"
print(sample)
h_cancer <- "H10"
h_tme <- "H3"
sce <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RDS_final/enhanced/", sample, "_sce_enhanced.rds"))
sce <- sce[,-1]
colData(sce)[, c(h_cancer, h_tme, "clusters")] <- df_all[df_all$sample==sample, c(h_cancer, h_tme, "clusters")]

v <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(1,2),], fill = h_cancer)
v2 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(4,5),], fill = h_tme)
v3 <- .make_triangle_subspots(colData(sce)[sce$clusters %in% c(3),], fill = "clusters")

ref_v <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/",sample,".rds"))
hires <- readRDS(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/",sample,".rds"))
for (spot in unique(v$spot)) {
  v$imagecol[v$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v$imagerow[v$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}


for (spot in unique(v2$spot)) {
  v2$imagecol[v2$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v2$imagerow[v2$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}

for (spot in unique(v3$spot)) {
  v3$imagecol[v3$spot == spot] <- ref_v$imagecol[ref_v$spot==spot] 
  v3$imagerow[v3$spot == spot] <- ref_v$imagerow[ref_v$spot==spot]
}

hires + geom_polygon(data=v,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn(h_cancer, colours = viridisLite::rocket(1000, alpha = 1, begin = 0, end = 1, direction = 1)[200:1000]) + 
  ggnewscale::new_scale("fill") + geom_polygon(data=v2,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~fill)) +  theme_void() + coord_equal()+
  scale_fill_gradientn(h_tme, colours = viridisLite::mako(1000, alpha = 1, begin = 0, end = 1, direction = 1)[200:1000]) + 
  ggnewscale::new_scale("fill") + geom_polygon(data=v3,  aes_(x=~imagerow+20, y=~imagecol+100, group=~spot, fill=~as.factor(fill))) +
  scale_fill_manual("Buffer", values = "gray41") 








ggsave(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/RF_plots/pos_RF/", sample, "_", h_cancer, "_", h_tme, ".png"), bg = "white", width = 7, height = 7)


    
