

STobject <- readRDS("Desktop/IJC/TFG/RDS/OV4A.rds")



df <- as.data.frame(t(STobject@assays$SCT@data))
df$id <- rownames(df)

STobject.images <- Read10X_Image("Desktop/IJC/datasets/IGTP/4A/spatial/", image.name = 'tissue_hires_image.png')
STobject.images@scale.factors$hires #0.9469697

rownames(STobject.images@coordinates) <- paste0(rownames(STobject.images@coordinates), "_1")
df2 <- STobject.images@coordinates[colnames(STobject), c("imagerow", "imagecol")]
colnames(df2) <- c("Y_hires_image", "X_hires_image")


df <- cbind(df, df2)
df <- df %>%
  select(id, everything())

write.table(df, "Desktop/df_input.tsv", sep = "\t", row.names = F, quote = F)


STobject <- readRDS("Desktop/IJC/TFG/RDS/OV4A.rds")



df <- as.data.frame(t(STobject@assays$SCT@data))
df$id <- rownames(df)

STobject.images <- Read10X_Image("Desktop/IJC/datasets/IGTP/4A/spatial/", image.name = 'tissue_hires_image.png')
STobject.images@scale.factors$hires #0.9469697

rownames(STobject.images@coordinates) <- paste0(rownames(STobject.images@coordinates), "_1")
df2 <- STobject.images@coordinates[colnames(STobject), c("imagerow", "imagecol")]
colnames(df2) <- c("Y_hires_image", "X_hires_image")


df <- cbind(df, df2)
df <- df %>%
  select(id, everything())

write.table(df, "Desktop/df_input_normalized.tsv", sep = "\t", row.names = F, quote = F)
