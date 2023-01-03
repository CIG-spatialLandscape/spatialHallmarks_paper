


x <- read.table("Desktop/IJC/datasets/IGTP/figuresPaper/tables/FunctionalEnrichment.xls", header = T, sep = "\t")
head(x)
colnames(x)[21] <- "Tumor type"
x <- x[,-20]
write.table(x, "Desktop/IJC/datasets/IGTP/figuresPaper/tables/Funct.xls", sep = "\t", row.names = F)


x <- read.table("Desktop/IJC/datasets/IGTP/figuresPaper/tables/hallmarks_stats.xls", header = T, sep = ",")


x$Hallmark <- sapply(x$Hallmark, function(h) {
  hallmark_names[as.numeric(substr(h, 2, nchar(h)))] 
})

write.table(x, "Desktop/IJC/datasets/IGTP/figuresPaper/tables/stats.xls", sep = "\t", row.names = F)
write.table(x, "Desktop/IJC/datasets/IGTP/figuresPaper/tables/tbl_S4.xls", sep = "\t", row.names = F)

