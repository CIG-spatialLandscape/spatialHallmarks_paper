STobject <- readRDS("Desktop/enhanced/new/C39_enhanced.rds")

VlnPlot(STobject, features = c("H3", "H4", "H6","H11", "H12", "H7"), group.by = "NMF_clusters")
p



embeddings <- STobject@reductions$NMF@cell.embeddings #extract factors embeddings
hallmarks <- STobject@meta.data[, paste0("H", 1:13)] #extract hallmarks scores
mat <- matrix(0, 11, 13) #create a matrix of 0s; dim (factors*hallmarks)
#name columns and rownames
colnames(mat) <- colnames(hallmarks)
rownames(mat) <- paste0(STobject@images[[1]]@key, "_", 1:11)
#fill each position of the matrix
for (h in 1:13) {
  for (factor in 1:11) {
    mat[factor, h] <- cor(embeddings[,factor], hallmarks[,h]) #compute the pearson correlation between factor activity and hallmarks
  }
}

mat <- as.data.frame(mat) #create a data frame from the matrix to add a column with string


pheatmap(t(mat), color =   rev(paletteer_c("ggthemes::Orange-Blue Diverging", 50)))



library(utils)

rforge <- "http://r-forge.r-project.org"

install.packages("estimate", repos=rforge, dependencies=TRUE)



library(estimate)
estimateScore()



#check how many signature genes are in each hallmark

stromal <- SI_geneset[1,2:142]
immune <- SI_geneset[2,2:142]
H <- read.table("Desktop/IJC/datasets/IGTP/figuresPaper/Hallmarks_genes_PthCommons2019_final_onlyGenes.txt", sep = "\t", header = T)
#Split genes by the corresponding hallmark 
H_features <- list("1"= H[H$H =="H1",1],"2"= H[H$H =="H2",1], "3"= H[H$H =="H3",1],"4"=H[H$H =="H4",1],"5"= H[H$H =="H5",1],"6"= H[H$H =="H6",1],"7"= H[H$H =="H7",1],
                   "8"= H[H$H =="H8",1],"9"= H[H$H =="H9",1],"10"= H[H$H =="H10",1],"11"= H[H$H =="H11",1],"12"= H[H$H =="H12",1], "13"= H[H$H =="H13",1])

table(stromal %in% H[[1]])
table(H_features$`1` %in% stromal)
table(H_features$`2` %in% stromal)
table(H_features$`3` %in% stromal)      
table(H_features$`4` %in% stromal)
table(H_features$`5` %in% stromal)
table(H_features$`6` %in% stromal)
table(H_features$`7` %in% stromal)
table(H_features$`8` %in% stromal)
table(H_features$`9` %in% stromal)
table(H_features$`10` %in% stromal)
table(H_features$`11` %in% stromal)
table(H_features$`12` %in% stromal)
table(H_features$`13` %in% stromal)



table(immune %in% H[[1]])
table(H_features$`1` %in% immune)
table(H_features$`2` %in% immune)
table(H_features$`3` %in% immune)      
table(H_features$`4` %in% immune)
table(H_features$`5` %in% immune)
table(H_features$`6` %in% immune)
table(H_features$`7` %in% immune)
table(H_features$`8` %in% immune)
table(H_features$`9` %in% immune)
table(H_features$`10` %in% immune)
table(H_features$`11` %in% immune)
table(H_features$`12` %in% immune)
table(H_features$`13` %in% immune)


#score stroma and immune using AddModuleScore

signature <- list(immune=t(immune), stromal=t(stromal))

STobject <- AddModuleScore(
  STobject,
  signature,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = NULL,
  name = "S",
  seed = 200,
  search = F
)

SpatialFeaturePlot(STobject, features=c("S1", "S2"), pt.size.factor=0.65)
VlnPlot(STobject, features = c("H3", "H6", "H7", "H1", "H11"), group.by = "NMF_clusters")


#estiamte package
library(phantasus)
expr <- ExpressionSet(as.matrix(GetAssayData(STobject, assay = "SCT")))
write.gct(expr, "Desktop/test.gct", gzip=F)

estimateScore(input.ds = "Desktop/test.txt", output.ds = "Desktop/out.gct", platform = "illumina")

a <- read.table("Desktop/out.gct", sep = "\t", skip = 3)


expr <- as.matrix(GetAssayData(STobject, assay = "SCT"))
write.table(expr, "Desktop/test.txt", quote = F, sep = "\t")

estimateScore(input.ds = "Desktop/test.txt", output.ds = "Desktop/out.gct", platform = "illumina")
a <- read.table("Desktop/out.gct", sep = "\t", skip = 3)
a$pseudo <- c(2200, 2400, 4600)
colnames(a) <- colnames(STobject)
STobject$StromalScore <- as.numeric(t(a[1,3:ncol(a)]))
STobject$ImmuneScore <- t(a[2,3:ncol(a)])
STobject$TumorScore <- t(a[3,3:ncol(a)])

SpatialFeaturePlot(STobject, features = c("StromalScore", "ImmuneScore", "TumorScore"), pt.size.factor = 0.65)
VlnPlot(STobject, features = c("StromalScore", "ImmuneScore", "TumorScore"), group.by = "NMF_clusters")
VlnPlot(STobject, features = "KRT19")

plotPurity("Desktop/out.gct")


#

pseudobulk <- AverageExpression(STobject, assays = "SCT", group.by = "NMF_clusters" )
pseudobulk <- pseudobulk[[1]]

#create GCT
df <- data.frame(NAME=rownames(pseudobulk), Description=rownames(pseudobulk))
df <- cbind(df, pseudobulk)
write.table(df, "Desktop/test39.gct", quote = F, sep = "\t", row.names = F)
estimateScore(input.ds = "Desktop/test39.gct", output.ds = "Desktop/out39.gct")
a <- read.table("Desktop/out39.gct", sep = "\t", skip = 3)
plot(1:11, as.numeric(t(a[1,3:ncol(a)])))
plot(1:11, as.numeric(t(a[2,3:ncol(a)])))
plot(1:11, as.numeric(t(a[3,3:ncol(a)])))
plot(1:11, as.numeric(t(a[4,3:ncol(a)])))
