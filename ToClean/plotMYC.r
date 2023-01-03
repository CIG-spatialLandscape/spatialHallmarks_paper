
#################### TF activity ###################
library(readr)
library(tidyr)
library(tibble)
library(stringr)
library(dplyr)
library(purrr)
library(viper)
library(progeny)
library(OmnipathR)

## Output: Object of class regulon. See viper package.
df2regulon = function(df) {
  regulon = df %>%
    split(.$tf) %>%
    map(function(dat) {
      tf = dat %>% distinct(tf) %>% pull()
      targets = setNames(dat$mor, dat$target)
      likelihood = dat$likelihood
      list(tfmode =targets, likelihood = likelihood)
    })
  return(regulon)
}

## <<DOROTHEA>>  

dorothea_regulon_human = read_csv("https://raw.githubusercontent.com/saezlab/ConservedFootprints/master/data/dorothea_benchmark/regulons/dorothea_regulon_human_v1.csv")

# We obtain the regulons based on interactions with confidence level A, B and C
DefaultAssay(STobject) <- "SCT"

regulon = dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C")) %>%
  dplyr::filter(tf == "MYC") %>%
  df2regulon()


tf_act_mat = viper(eset = as.matrix(STobject[["SCT"]]@data),
                   regulon = regulon, nes = TRUE,
                   method = "scale", minsize = 4,
                   eset.filter = FALSE,
                   verbose = FALSE, cores = 14)  
STobject$MYC <- NULL
STobject$MYC <- tf_act_mat

## Repeated regulon RFXAP
#tf_act_mat = tf_act_mat[!duplicated(tf_act_mat),]
STobject[['dorothea']] <- NULL
STobject[['dorothea']] = CreateAssayObject(counts = tf_act_mat)
DefaultAssay(STobject) <- "dorothea"
SpatialFeaturePlot(STobject, features = "MYC")
rm(tf_act_mat)



p <- SpatialFeaturePlot(STobject, features = "MYC", pt.size.factor = 0.65) + scale_fill_gradientn("", colours = viridis::inferno(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 40)) + ggtitle("MYC")

pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/MYC.pdf")
plot(p)
dev.off()


model <- progeny::getModel(organism = "Human", top = 1000)
common_genes <- intersect(rownames(GetAssayData(STobject, assay = "SCT")), rownames(model))
progeny_scores <- t(model)[, common_genes] %*% GetAssayData(STobject, assay = "SCT")[common_genes, ]
progeny_scores <- scale(t(as.matrix(progeny_scores)))

STobject[['progeny']] <- CreateAssayObject(counts = t(progeny_scores))
DefaultAssay(STobject) <- "progeny"

p <- SpatialFeaturePlot(STobject, features = "TGFb", pt.size.factor = 0.65) + scale_fill_gradientn("", colours = viridis::inferno(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 40)) + ggtitle("TGFb") + scale_color_manual("red")

pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/TGFb.pdf")
plot(p)
dev.off()

###########3

palette <- RColorBrewer::brewer.pal(12, name = "Paired")
palette <- RColorBrewer::brewer.pal(12, name = "Paired")[c(6, 8, 11, 1, 4)]
