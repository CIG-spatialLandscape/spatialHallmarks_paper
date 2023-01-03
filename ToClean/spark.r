library(SPARK)
setwd("FFPE_NP")
p = 'Visium_FFPE_Human_Normal_Prostate_filtered_feature_bc_matrix.h5'

expr.data <- Seurat::Read10X_h5(p)

info <- cbind.data.frame(x=as.numeric(sapply(strsplit(colnames(expr.data), split="x"), "[", 1)),
                         y=as.numeric(sapply(strsplit(colnames(expr.data), split="x"), "[", 2)),
                         total_counts=apply(expr.data,2,sum))

img <- Seurat::Read10X_Image("Visium_FFPE_Human_Normal_Prostate_spatial/spatial")
Seurat::DefaultAssay(img) <-"Spatial"

info$x <- img@coordinates[4]
info$y <- img@coordinates[5]

spark <- CreateSPARKObject(counts=expr.data, location=info[,1:2], percentage=0.1, min_total_counts=10)

spark@lib_size <- apply(spark@counts, 2, sum)

spark <- spark.vc(spark, covariates=NULL, lib_size=spark@lib_size, num_core=1, verbose=T)


spark <- spark.test(spark, 
                    check_positive = T, 
                    verbose = T)

head(spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")])
