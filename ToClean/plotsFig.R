hires <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/image/Colorectal.rds")
v <- readRDS("Desktop/IJC/datasets/IGTP/figuresPaper/hiresplot/bayes/Colorectal.rds")
v <- v[v$spot != "subspot_1.1",]
df_all$spotid <- sapply(rownames(df_all), function(spot){
  x <- strsplit(spot,  ".", fixed = T)
  paste0(x[[1]][1], ".", substr(x[[1]][2],1, 1))
})
sub_data <- df_all[df_all$sample == "Colorectal",]
sub_data <- read.table("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/output_df/Colorectal.txt")

library(paletteer)
d <- data.frame(estimate=sub_data$estimate, score=sub_data$score)
model <- lm(estimate~score, data = d)
sub_data$residuals <- model$residuals
library(paletteer)

cluster <- c("5" = "lightgoldenrod1", "4" = "lightgoldenrod3", "3" = "lightpink2", "2"= "orchid3", "1"= "orchid4")


v$fill <- sub_data[v$spot,"estimate"]
hires + geom_polygon(data=v,  aes_(x=~imagerow, y=~imagecol, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  #scale_fill_viridis_c()
  #scale_fill_gradient2(low = "orchid4", mid = "lightpink2", high = "lightgoldenrod1")
  scale_fill_gradientn(colours = rev(cluster))

v$fill <- sub_data[v$spot,"clusters"]
hires + geom_polygon(data=v,  aes_(x=~imagecol, y=~imagerow, group=~spot, fill=~factor(fill))) +  theme_void() + coord_equal() +
  scale_fill_manual(values = rev(c("lightgoldenrod1", "lightgoldenrod3", "lightpink2",  "orchid3",  "orchid4"))) + labs(fill="") 

v$fill <- sub_data[v$spot,"residuals"]
hires + geom_polygon(data=v,  aes_(x=~imagecol, y=~imagerow, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  #scale_fill_viridis_c()
  scale_fill_gradient2(low = "navy", mid = "gainsboro", high = "firebrick3")
  #scale_fill_gradientn(colours = paletteer_c("grDevices::Purple-Green", 30) )

ggplot(sub_data, aes(x=score, y=estimate, col=residuals)) + geom_point() + geom_smooth(method = "lm", color="black") +
scale_color_gradient2(low = "navy", mid = "gainsboro", high = "firebrick3") + theme_classic() + labs(x="NEIGHBORHOOD score", y="ESTIMATE score", col="Residuals") + 
   scale_x_continuous(labels = c("Cancer", "TME"), breaks = c(-50, 100)) + 
   scale_y_continuous(labels = c("Cancer", "TME"), breaks = c(-5000, 10000)) + 
   theme(axis.ticks = element_blank(), axis.text = element_text(size = 13), axis.title = element_text(size = 17)) 


df <- read.table("Desktop/df_tme58.txt")

df$spotid <- sapply(rownames(df), function(spot){
  x <- strsplit(spot,  ".", fixed = T)
  paste0(x[[1]][1], ".", substr(x[[1]][2],1, 1))
})

tmp_v <- v
tmp_v <-  tmp_v[tmp_v$spot %in% df$spotid[df$sample=="Colorectal"],]
tmp_v$fill <- df[df$sample=="Colorectal",][tmp_v$spot,"H11_cancer"]
hires + geom_polygon(data=tmp_v,  aes_(x=~imagecol, y=~imagerow, group=~spot, fill=~fill)) +  theme_void() + coord_equal() +
  #scale_fill_viridis_c()
  scale_fill_gradient2(low = "navy", mid = "gainsboro", high = "firebrick3")
#scale_fill_gradientn(colours = paletteer_c("grDevices::Purple-Green", 30) )

###
lm_mat.basal$cnv <- "Monoclonal"
lm_mat.basal$cnv[lm_mat.basal$Sample %in% lowdiff_samples] <- "Low"
lm_mat.basal$cnv[lm_mat.basal$Sample %in% highdiff_samples] <- "High"

ggplot(lm_mat.basal, aes(x=cnv, fill=cnv, y=RsquaredCancer)) + geom_boxplot()
lm_mat.basal$y <- lm_mat.basal$RsquaredCancer + lm_mat.basal$diffCancer
ggboxplot(lm_mat.basal, x="cnv", fill="cnv", y="y") + stat_compare_means(comparisons = list(c("High", "Low"), c("Monoclonal", "Low"), c("High", "Monoclonal"))) + 
  labs(y="Rsquared")

ggboxplot(lm_mat.basal, x="cnv", fill="cnv", y="RsquaredCancer") + stat_compare_means(comparisons = list(c("High", "Low"), c("Monoclonal", "Low"), c("High", "Monoclonal")))
ggboxplot(lm_mat.basal, x="cnv", fill="cnv", y="RsquaredTME") + stat_compare_means(comparisons = list(c("High", "Low"), c("Monoclonal", "Low"), c("High", "Monoclonal")))
ggboxplot(lm_mat.basal, x="cnv", fill="cnv", y="RsquaredBuffer") + stat_compare_means(comparisons = list(c("High", "Low"), c("Monoclonal", "Low"), c("High", "Monoclonal")))
ggboxplot(lm_mat.basal, x="cnv", fill="cnv", y="Rsquared") + stat_compare_means(comparisons = list(c("High", "Low"), c("Monoclonal", "Low"), c("High", "Monoclonal")))


df_hallmarks <- df_all[, c(paste0("H", 1:13),"clusters")]
colnames(df_hallmarks)[1:13] <- hallmark_names
df_hallmarks.long <- melt(df_hallmarks, id.vars = "clusters")

# Create Boxplot with a line plot using mean values
ggplot(df_hallmarks.long, aes(x=as.factor(clusters), y=value, fill=as.factor(clusters))) + 
  geom_boxplot(outlier.size = 0.2) + facet_wrap(~variable, scales = "free") + theme_classic2() + 
  scale_fill_manual(values=rev(c("lightgoldenrod1", "lightgoldenrod3", "lightpink2",  "orchid3",  "orchid4"))) + 
  labs(fill="ESTIMATE clusters", y="Hallmark activity", x="")


#####
