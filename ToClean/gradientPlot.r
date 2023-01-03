library(pheatmap)
library(paletteer)
library(grid)

sub <- fitting[fitting$sample=="Ovarian" & fitting$hallmark == "H11" & fitting$compartment=="Cancer",]
sub1 <- fitting[fitting$sample=="OV4A" & fitting$hallmark == "H11" & fitting$compartment=="Cancer",]

min(sub$x)

#plot 1 sample different hallmarks
sub <- fitting[fitting$sample=="OV4A" & fitting$compartment=="Cancer",]


sub <- reshape(sub[c("x", "y", "hallmark")], idvar = "x",timevar = "hallmark", direction = "wide")

mat <- data.frame(sub, row.names = 1)


pheatmap(t(mat),
         cluster_cols = F, cluster_rows = T,
         scale = "row",
         show_colnames = F,
         color =   rev(paletteer_c("ggthemes::Orange-Blue Diverging", 50) ))

#plot 2 sample 1 hallmark
sub <- fitting[fitting$sample %in% c("Ovarian", "OV4A") & fitting$hallmark == "H11" & fitting$compartment=="Cancer",]
sub$x <- round(sub$x)
sub$x <- as.character(sub$x)
sub <- reshape(sub[c("x", "y", "sample")], idvar = "x", timevar = "sample", direction = "wide")

mat <- data.frame(sub, row.names = 1)
mat <- mat[as.character(sort(as.numeric(rownames(mat)))),]

pheatmap(t(mat),
         cluster_cols = F, cluster_rows = F,
         scale = "row",
         show_colnames = F,
         color =   rev(paletteer_c("ggthemes::Orange-Blue Diverging", 50) ), na_col = "black")



#plot 2 sample 1 hallmark

sub <- fitting[!fitting$sample %in% c("Glioblastoma", "Cervical", "Acinar", "IC") & fitting$hallmark == "H11" & fitting$compartment=="Cancer",]
sub <- fitting2[fitting2$hallmark == "H12",]


sub$x <- round(sub$x)
maximum <- paste0(max(sub$x), "μm")
minimum <- paste0(min(sub$x), "μm")

sub$x <- as.character(sub$x)
sub <- reshape(sub[c("x", "y", "sample")], idvar = "x", timevar = "sample", direction = "wide")

mat <- data.frame(sub, row.names = 1)
mat <- mat[as.character(sort(as.numeric(rownames(mat)))),]


ann <- data.frame(row.names = colnames(mat))
ann$Site <- sapply(rownames(ann), function(sample){
  annotation_tumor[[substr(sample, 3, nchar(sample))]]
})


pdf("Desktop/test.pdf")
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=1, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(t(mat),
         cluster_cols = F, cluster_rows = T,
         scale = "none",
         show_colnames = F,
         color =   rev(paletteer_c("ggthemes::Orange-Blue Diverging", 50) ), na_col = "black")
setHook("grid.newpage", NULL, "replace")
grid.text(minimum, y=-0.02, x=0.15,  gp=gpar(fontsize=12))
grid.text(maximum, y=-0.02, x=0.6,  gp=gpar(fontsize=12))
dev.off()

write.table(fitting2, "Desktop/distance/fitted/matrix_all.txt", sep = "\t")
fitting <- fitting2
for (compartment in unique(fitting$compartment)) {
  for (hallmark in unique(fitting$hallmark)) {
    sub <- fitting[fitting$hallmark == "H2",]
    #sub$x <- round(sub$x)
    #maximum <- paste0(max(sub$x), "μm")
    #minimum <- paste0(min(sub$x), "μm")
    #sub$x <- as.character(sub$x)
    sub <- reshape(sub[c("x", "y", "sample")], idvar = "x", timevar = "sample", direction = "wide")
    mat <- data.frame(sub, row.names = 1)
    mat <- mat[as.character(sort(as.numeric(rownames(mat)))),]
    ann <- data.frame(row.names = colnames(mat))
    ann$Site <- sapply(rownames(ann), function(sample){
      annotation_site[[substr(sample, 3, nchar(sample))]]
    })
    pdf(paste0("Desktop/heatmap/", compartment, "_", hallmark, ".pdf"))
    setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=1, height=0.9, name="vp", just=c("right","top"))), action="prepend")
    pheatmap(t(mat),
             cluster_cols = F, cluster_rows = T, annotation_row = ann,
             scale = "row",
             show_colnames = T,
             color =   rev(paletteer_c("ggthemes::Orange-Blue Diverging", 50) ), na_col = "black")
    setHook("grid.newpage", NULL, "replace")
    grid.text(minimum, y=-0.02, x=0.15,  gp=gpar(fontsize=12))
    grid.text(maximum, y=-0.02, x=0.6,  gp=gpar(fontsize=12))
    dev.off()
  }
}
