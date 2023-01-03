summary_comparisions <- read.table("Desktop/df_proximity.txt", header = T, sep = "\t")

summary_comparisions$Cancer[summary_comparisions$Cancer=="Breast"] <- "Invasive Lobular Carcinoma (Breast)"
summary_comparisions$Cancer[summary_comparisions$Cancer=="Ductal"] <- "Invasive Ductal Carcinoma (Breast)"
summary_comparisions$Cancer[summary_comparisions$Cancer=="IC"] <- "Adenocarcinoma (Prostate)"
summary_comparisions$Cancer[summary_comparisions$Cancer=="OV4A"] <- "HGSOC (Ovary)"
summary_comparisions$Cancer[summary_comparisions$Cancer=="Ovarian"] <- "Endometrial Adenocarcinoma (Ovary)"
summary_comparisions$Cancer[summary_comparisions$Cancer=="Acinar"] <- "Acinar Cell Carcinoma (Prostate)"
summary_comparisions$Cancer[summary_comparisions$Cancer=="Colorectal"] <- "Invasive Adenocarcinoma (Colorectal)"

summary_comparisions$unique <- paste0(summary_comparisions$Cancer," ", summary_comparisions$Cancer_factor)

summary_comparisions$unique[summary_comparisions$unique=="Invasive Lobular Carcinoma (Breast) 4"] <- "Invasive Lobular Carcinoma (Breast) - Region 1"
summary_comparisions$unique[summary_comparisions$unique=="Invasive Lobular Carcinoma (Breast) 5"] <- "Invasive Lobular Carcinoma (Breast) - Region 2"
summary_comparisions$unique[summary_comparisions$unique=="Invasive Ductal Carcinoma (Breast) 1"] <- "Invasive Ductal Carcinoma (Breast) - Region 1"
summary_comparisions$unique[summary_comparisions$unique=="Invasive Ductal Carcinoma (Breast) 2"] <- "Invasive Ductal Carcinoma (Breast) - Region 2"
summary_comparisions$unique[summary_comparisions$unique=="Adenocarcinoma (Prostate) 1"] <- "Adenocarcinoma (Prostate) - Region 1"
summary_comparisions$unique[summary_comparisions$unique=="Adenocarcinoma (Prostate) 5"] <- "Adenocarcinoma (Prostate) - Region 2"
summary_comparisions$unique[summary_comparisions$unique=="Invasive Adenocarcinoma (Colorectal) 1"] <- "Invasive Adenocarcinoma (Colorectal) - Region 1"
summary_comparisions$unique[summary_comparisions$unique=="Invasive Adenocarcinoma (Colorectal) 4"] <- "Invasive Adenocarcinoma (Colorectal) - Region 2"
summary_comparisions$unique[summary_comparisions$unique=="HGSOC (Ovary) 1"] <- "HGSOC (Ovary)"
summary_comparisions$unique[summary_comparisions$unique=="Endometrial Adenocarcinoma (Ovary) 1"] <-  "Endometrial Adenocarcinoma (Ovary)"
summary_comparisions$unique[summary_comparisions$unique=="Acinar Cell Carcinoma (Prostate) 2"] <- "Acinar Cell Carcinoma (Prostate) 2"
summary_comparisions$unique[summary_comparisions$unique=="Glioblastoma 3"] <- "Glioblastoma"

summary_comparisions$Diff <- as.numeric(summary_comparisions$Diff)
summary_comparisions$Sgnif <- as.numeric(summary_comparisions$Sgnif)


summary_comparisions$Diff[summary_comparisions$Sgnif > 0.0001] <- 0
summary_comparisions$Diff[summary_comparisions$Diff < 0.2 & summary_comparisions$Diff > -0.2] <- 0
summary_comparisions$color <- "a"

summary_comparisions$color[summary_comparisions$Diff >= 0 & summary_comparisions$Hallmark %in% hallmark_names[c(2,4,8,9,11,12)]] <- "TME-proximal"
summary_comparisions$color[summary_comparisions$Diff < 0 & summary_comparisions$Hallmark %in% hallmark_names[c(2,4,8,9,11,12)]] <- "TME-distant"

summary_comparisions$color[summary_comparisions$Diff >= 0 & summary_comparisions$Hallmark %in% hallmark_names[c(1,3,5,6,7,10,13)]] <- "Cancer-proximal"
summary_comparisions$color[summary_comparisions$Diff < 0 & summary_comparisions$Hallmark %in% hallmark_names[c(1,3,5,6,7,10,13)]] <- "Cancer-distant"



#summary_comparisions$unique <- factor(summary_comparisions$unique, levels=colnames(mat)[c$order])


summary_comparisions$group_hallmark <- "TME Hallmarks"
summary_comparisions$group_hallmark[summary_comparisions$Hallmark %in% hallmark_names[c(2,4,8,9,11,12)]] <- "Cancer Hallmarks"
summary_comparisions$group_hallmark <- factor(summary_comparisions$group_hallmark, levels=c("TME Hallmarks", "Cancer Hallmarks"))


                
palette1 <- RColorBrewer::brewer.pal(9,name="YlOrBr")
palette2 <- RColorBrewer::brewer.pal(9 , name = "YlGnBu")

l <- names(sort(table(summary_comparisions$color[summary_comparisions$Diff != 0], summary_comparisions$unique[summary_comparisions$Diff != 0])["TME-proximal",], decreasing = F))
summary_comparisions$unique <- factor(summary_comparisions$unique, levels=l)
l2 <- names(sort(table(summary_comparisions$colo[summary_comparisions$Diff != 0], summary_comparisions$Hallmark[summary_comparisions$Diff != 0])[c("TME-proximal"),hallmark_names[c(2,4,8,9,11,12)]], decreasing = T))
l3 <- names(sort(table(summary_comparisions$color[summary_comparisions$Diff != 0], summary_comparisions$Hallmark[summary_comparisions$Diff != 0])[c("Cancer-proximal"),hallmark_names[c(1,3,5,6,7,10,13)]],decreasing = F))

summary_comparisions$Hallmark <- factor(summary_comparisions$Hallmark, levels=c(l2,l3))



p <- ggplot(summary_comparisions[summary_comparisions$Diff != 0, ], 
       aes(x=Hallmark, y=as.factor(unique), size=abs(Diff), color=as.factor(color))) + geom_point() +
  theme_classic() + theme(axis.text.x = element_text(angle = -45, hjust = 0)) + scale_color_manual(values = c("#F46D43", "#A50026","#74ADD1", "#313695")) +
  facet_grid(cols=vars(group_hallmark), drop = T, scales = "free") + labs(size="Mean difference in proximity", color = "Proximity", y = "Cancer regions across cancer types") + scale_size(breaks = c(0.01, 0.3, 1.4))
p
pdf("Desktop/IJC/datasets/IGTP/figuresPaper/figures/correlationheatmap.pdf", width = 12, height = 5)
plot(p)
dev.off()    



tidyr::spread(data = summary_comparisions, key = c("unique", "Hallmark"), value = "Diff")

mat <- acast(summary_comparisions[,c("Hallmark", "unique", "Diff")], unique ~ Hallmark)
mat[ mat > 0.3 ] <- 1
mat[ mat < 0.3 & mat > -0.3] <- 0
mat[ mat < -0.3 ] <- -1

min(mat)
pheatmap(mat, cluster_cols = T)

mat <- apply(mat, 1, scale)
colnames(mat) <- NULL
d <- dist(ma, ...))
names(d) <- colnames(mat)
rownames(d) <- rownames(mat)

c <- hclust(dist(t(mat)))
c$labels <- colnames(mat)
plot(c)

c$order
colnames(mat)[c$order]
