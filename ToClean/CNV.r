

STobject <- readRDS("Desktop/IJC/datasets/Public/HCC-1T/RDS/HCC-1T.rds")
STobject$CNV <- labels[names,2]


mat <- as.matrix(STobject@assays$RNA@counts)
write.table(t(mat), "Desktop/CNV/mat.txt")
  

labels <- read.csv("Desktop/CNV/out/labels_name.csv")
states <- read.csv("Desktop/CNV/out/states_name.csv", row.names = 1)



names <- paste0(STobject@images[[1]]@coordinates[,"row"], ".0x",STobject@images[[1]]@coordinates[,"col"], ".0")


table(names %in% labels$X)
rownames(labels) <- labels$X


labels[names,2]





library(ggplot2)

ggplot(states, aes(x=X, y=X1)) +  geom_jitter()



pheatmap(t(states), cluster_cols = F, cluster_rows = F, show_colnames = F)




gtf <- rtracklayer::import('Desktop/IJC/datasets/genes.gtf')
gtf_df=as.data.frame(gtf)
