##################################################
## Project: Cancer Hallmarks
## Script purpose: Use BayesSpace enhancing algorithm in Visium Spatial transcriptomics data
## Author: Sergi Cervilla* & Mustafa Sibai*
##################################################
library(BayesSpace)

#signle cell experiment object at spot level
input_file <- "" 
output_file <- ""
ST_sce <- readRDS(input_file)  
#check image orientation - swap imagerow and imagecol
if (cor(ST_sce$imagecol, ST_sce$col) < 0.9) {
  tmp <- ST_sce$imagecol
  ST_sce$imagecol <- ST_sce$imagerow
  ST_sce$imagerow <- tmp
}
#enhance single cell experiment object
ST_sce.enhanced <- spatialEnhance(ST_sce, q=length(unique(ST_sce$spatial.cluster)), d=15, platform="Visium", init=ST_sce$spatial.cluster, 
                                  nrep=200000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE)
#save enhanced object (without gene imputation)
saveRDS(ST_sce.enhanced, output_file)
