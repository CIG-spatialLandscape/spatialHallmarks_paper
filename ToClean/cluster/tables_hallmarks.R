genes <- read.table("Desktop/IJC/datasets/IGTP/figuresPaper/Genes/genes_tmp_input_Acinar.gct.txt", header = F)

hallmark_genes <- read.table("Desktop/IJC/HAGs_intersect.txt", header = T)
#filter genes that are present in estimate
estimate_genes <- c(as.character(estimate::SI_geneset[1,-1]), as.character(estimate::SI_geneset[2,-1]))
hallmark_genes <- hallmark_genes[!hallmark_genes$gene %in% estimate_genes, ]

#genes detected in each case
genes <- read.table("Desktop/IJC/datasets/IGTP/figuresPaper/Genes/genes_tmp_input_Acinar.gct.txt", header = F)


samples <- list.files("Desktop/IJC/datasets/IGTP/figuresPaper/SCD/")
samples <- sapply(samples, function (x) {strsplit(x, split = ".", fixed = T)[[1]][1]})
df_hallmarks <- data.frame()
for (sample in samples) {
  genes <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/Genes/genes_tmp_input_", sample, ".gct.txt"), header = F)
  for (hallmark in paste0("H", 1:13)) {
   tbl <- table(hallmark_genes$gene[hallmark_genes$H == hallmark] %in% genes$V1) 
   genes_found <- hallmark_genes$gene[hallmark_genes$H == hallmark][hallmark_genes$gene[hallmark_genes$H == hallmark] %in% genes$V1]
   df_hallmarks <- rbind(df_hallmarks, c(sample, 
                     hallmark_names[[hallmark]],
                     as.numeric(tbl["TRUE"]),
                     as.numeric(tbl["FALSE"]),
                     round(as.numeric(tbl["TRUE"]/sum(tbl)),4)*100,
                     paste0(genes_found, collapse = ";")))
  }  
}
colnames(df_hallmarks) <- c("Sample", "Hallmark", "Genes Present", "Genes not Present", "Percentage", "Genes found")

df_hallmarks$Tissue <- sapply(df_hallmarks$Sample, function(x){tissue[[x]]})

#Plot
color_codes <-  list("Sustaining Proliferative Signaling"="#15CE59", 
                     "Evading Growth Suppressors" = "#701717",
                     "Avoiding Immune Destruction" = "#CB3BBD",
                     "Enabling Replicative Immortality" = "#4969D4",
                     "Tumor-Promoting Inflammation" = "#E6880D",
                     "Activating Invasion and Metastasis" = "#000000",
                     "Inducing Angiogenesis" = "#EE0F16",
                     "Genome Instability and Mutation" = "#132892",
                     "Resisting Cell Death" = "#8E909B",
                     "Deregulating Cellular Energetics" = "#71189E",
                     "Senescent cells" = "#05F3EB",
                     "Nonmutational Epigenetic reprogramming" = "#890269",
                     "Unlocking Phenotypic Plasticity" = "#95641A")
palette <- do.call(rbind, color_codes)
df_hallmarks$Hallmark <- factor(df_hallmarks$Hallmark, levels = names(color_codes))
ggplot(df_hallmarks, aes(x=Hallmark, y=as.numeric(Percentage), fill=Hallmark)) + 
  geom_boxplot() + ylim(c(0, 100)) + scale_fill_manual(values = palette) + theme_classic() + 
  labs(x="", y="Percentage of hallmark genes found") + theme(axis.text.x = element_text(angle = -45, hjust = 0),
                                                             axis.text = element_text(size = 15), 
                                                             axis.title.y = element_text(size=20))
ggplot(df_hallmarks, aes(x=Tissue, y=as.numeric(Percentage), fill=Tissue)) + 
  geom_boxplot() + ylim(c(0, 100)) + theme_classic() + 
  labs(x="", y="Percentage of hallmark genes found") + theme(axis.text.x = element_text(angle = -45, hjust = 0),
                                                             axis.text = element_text(size = 15), 
                                                             axis.title.y = element_text(size=20))

#
df_overlap <- data.frame()
for (hallmark in paste0("H", 1:13)) {
  l_genes <- lapply(df_hallmarks[df_hallmarks$Hallmark == hallmark_names[[hallmark]], "Genes found"], function(x) {str_split(x, pattern = ";")[[1]]})
  
  df_overlap <- rbind(df_overlap, c(hallmark_names[[hallmark]], round(length(Reduce(intersect, l_genes)) /length(hallmark_genes[hallmark_genes$H==hallmark, "gene"]), 4)))
}
colnames(df_overlap) <- c("Hallmark", "Overlap")
df_overlap$Hallmark <- factor(df_overlap$Hallmark, levels = names(color_codes))
ggplot(df_overlap, aes(x=Hallmark, y=as.numeric(Overlap)*100, fill=Hallmark)) + geom_bar(stat = "identity") +
  ylim(c(0, 100)) + scale_fill_manual(values = palette) + theme_classic() + 
  labs(x="", y="Percentage of overlap hallmark genes") + theme(axis.text.x = element_text(angle = -45, hjust = 0),
                                                             axis.text = element_text(size = 13), 
                                                             axis.title.y = element_text(size=16))
