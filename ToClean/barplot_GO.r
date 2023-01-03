library(dplyr)
library(ggpubr)
library(colorspace)
library(RColorBrewer)
library(gprofiler2)


top_30 <- markers %>% top_n(30, avg_log2FC)

top.genes <- list()
nclusters <- 5
for (cluster in 1:nclusters) {
  top.genes[[as.character(cluster)]] <- top_30[top_30$cluster==cluster,]$gene
}


pathways <- lapply(seq_along(top.genes), function(i) {
  gset <- top.genes[[i]]
  df <- NULL
  if (length(top.genes[[i]]) !=0) {
    df <- gost(query = gset, organism = "hsapiens", sources = "GO:BP", evcodes = TRUE, )$result
  }
  #df2 <- gost(query = gset, organism = 'gp__aJd5_EBf2_fFc', evcodes = TRUE)$result     #hallmarks from DBsig
  #df <- rbind(df1, df2)
  if (is.null(df)) return(NULL)
  df$cluster <- paste0("cluster_", i)
  return(df)
})


select_top <- function(st, factor) {
  vals <- st@reductions$NMF@feature.loadings[, factor]
  logvals <- log(vals)
  logvals <- logvals[is.finite(logvals)]
  thr <- mean(logvals) + 1.4*sd(logvals)
  names(logvals[logvals > thr])
}

loadings <- select_top(STobject, 3)

f3 <- gost(query = loadings, organism = "hsapiens", sources = "GO:BP", evcodes = TRUE, )$result
f3$cluster <- "cluster_3"




pathways <- do.call(rbind, pathways)
pathways <- rbind(pathways, f3)

pathways$GeneRatio <- pathways$intersection_size/pathways$query_size
#select top 30 (by significance)
pathways_top5 <- pathways %>% 
  group_by(cluster) %>%
  top_n(n = 3, wt = -log10(p_value))

pathways_top5 <- pathways_top5[-3,]





pathways_top5$cluster <- gsub("cluster_", "", pathways_top5$cluster)
# barplot for top5 BPs
palette <- RColorBrewer::brewer.pal(12, name = "Paired")[c(6, 8, 11, 1, 4)]
palette
#BPs_top20 <- arrange(BPs_top20, desc(cluster))

x <- ggbarplot(pathways_top5, x = "term_name", y = "GeneRatio",
               color = "cluster",  fill = "cluster",
               palette = palette,
               label = FALSE, lab.pos = "in", lab.col = "white",
               ggtheme = theme_pubclean(),
               #x.text.angle = 70,
               sort.by.groups = F,
               rotate = F, xlab = "GO:BP term"
) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  font("y.text", size = 11) +
  font("title", color = darken("#00AFBB", amount = 0.3)) +
  font("caption", face = "italic") +
  font("legend.text", size = 13) + scale_fill_manual(labels = c("Cancer", "Immune", "ECM", "ECM (Endothelial)", "ECM"), values = palette) + 
  scale_color_manual(labels = c("Cancer", "Immune", "ECM", "ECM (Endothelial)", "ECM"), values = palette) + 
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = -45, hjust = 0, size = 15), plot.margin = margin(0.5,6.5,0.5,0.5, "cm"))



x


