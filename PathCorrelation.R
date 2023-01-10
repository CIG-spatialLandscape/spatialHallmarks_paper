##################################################
## Project: Cancer Hallmarks
## Script purpose: Correlate hallmark activity and pathway scores across all samples
## Date: 22/12/2022
## Author: Sergi Cervilla & Mustafa Sibai
##################################################

library(dplyr)
library(ggplot2)
library(ggpubr)
library(colorspace)

### list files with pathway scores
files <- list.files("", pattern = ".txt" ,recursive = F)
#Pan-Cancer pathways within TME sub-spots
Pan_TME <- data.frame()
#Pan-Cancer pathways within Cancer sub-spots
Pan_cancer <- data.frame()

#for loop to load the pathway scores for each sample
for (file in files) {
  print(file)
  #load pathway score for a given sample
  sample <- read.delim("", header = T, check.names = F)
  #sample name
  sample.col <- data.frame(sample = rep(strsplit(file, split = ".txt")[[1]][1], nrow(sample)))
  #data frame containing sample name and pathway scores 
  sample <- cbind(sample.col, sample)
  #add tumor anatomic site metadata
  if (unique(sample.col$sample) %in% c("DU2", "DU3", "DU8", "DU12", "DU13")) {
    tumor.col <- data.frame(tumor = rep("Bladder", nrow(sample)))
    sample <- cbind(tumor.col, sample)
  } else if (unique(sample.col$sample) %in% c("Glioblastoma", "UKF242T", "UKF260T", "UKF269T", "UKF275T")) {
    tumor.col <- data.frame(tumor = rep("Brain", nrow(sample)))
    sample <- cbind(tumor.col, sample)
  } else if (unique(sample.col$sample) %in% c("Breast", "BreastA", "Ductal", "DuctalFFPE", "TNBCA", "M1", "M2", "M3", "M4")) {
    tumor.col <- data.frame(tumor = rep("Breast", nrow(sample)))
    sample <- cbind(tumor.col, sample)
  } else if (unique(sample.col$sample) %in% c("cHC1T", "HCC1T", "HCC2T", "HCC5D", "ICC1L")) {
    tumor.col <- data.frame(tumor = rep("Liver", nrow(sample)))
    sample <- cbind(tumor.col, sample)
  } else if (unique(sample.col$sample) %in% c("P259_H2A2", "P264", "P270", "P288", "P306")) {
    tumor.col <- data.frame(tumor = rep("Pancreas", nrow(sample)))
    sample <- cbind(tumor.col, sample)
  } else if (unique(sample.col$sample) %in% c("OV4A", "OVD1", "CUP295", "OVFFPE", "Ovarian")) {
    tumor.col <- data.frame(tumor = rep("Ovary", nrow(sample)))
    sample <- cbind(tumor.col, sample)
  } else if (unique(sample.col$sample) %in% c("Colorectal", "CRC1", "CRC2", "Intestine", "Co1", "Co2", "Co3", "Co4")) {
    tumor.col <- data.frame(tumor = rep("Colorectal", nrow(sample)))
    sample <- cbind(tumor.col, sample)
  } else if (unique(sample.col$sample) %in% c("Acinar", "IC", "PC1", "PC2")) {
    tumor.col <- data.frame(tumor = rep("Prostate", nrow(sample)))
    sample <- cbind(tumor.col, sample)
  } else if (unique(sample.col$sample) %in% c("P1", "P3", "P4", "P5", "P6", "P7", "P8")) {
    tumor.col <- data.frame(tumor = rep("Lung", nrow(sample)))
    sample <- cbind(tumor.col, sample)
  } else {
    tumor.col <- data.frame(tumor = rep("Kidney", nrow(sample)))
    sample <- cbind(tumor.col, sample)
  }
  #change rownames ID
  rownames(sample) <- 1:nrow(sample)
  
  #Subset TME sub-spots (include Buffer region)
  sample_TME <- filter(sample, estimate.cluster %in% c(3,4,5))
  #Scale pathway scores within the subset
  sample_TME[,-c(1:3)] <- scale(sample_TME[,-c(1:3)])
  #add TME sub-spots to the rest of samples
  Pan_TME <- dplyr::bind_rows(Pan_TME, sample_TME)
  #Subset Cancer sub-spots (include Buffer region)
  sample_cancer <- filter(sample, estimate.cluster %in% c(1,2,3))
  #Scale pathway scores within the subset
  sample_cancer[,-c(1:3)] <- scale(sample_cancer[,-c(1:3)])
  #add Cancer sub-spots to the rest of samples
  Pan_cancer <- dplyr::bind_rows(Pan_cancer, sample_cancer)
  #clean memory
  gc()
}

#clean memory
rm(sample, sample_cancer, sample_TME, sample.col, tumor.col)
gc()

## TME
#create a list of empty data frames for each hallmark
corr.TME <- list(H1 = data.frame(), H2= data.frame(), H3 = data.frame(),
                 H4 = data.frame(), H5 = data.frame(), H6 = data.frame(),
                 H7 = data.frame(), H8 = data.frame(), H9 = data.frame(),
                 H10 = data.frame(), H11 = data.frame(), H12 = data.frame(), H13 = data.frame())

#for each hallmark add the correlation (coefficient and significance) of hallmark activity and pathway score 
for (hallmark in paste0("H", 1:13)) {
  tmp <- sapply(1:ncol(Pan_TME[,c(hallmark, grep(paste0(hallmark, "_"), colnames(Pan_TME[,-c(1:16)]), value = T))][,-1]), function(x) {
    data.frame(pathway = colnames(Pan_TME[,c(hallmark, grep(paste0(hallmark, "_"), colnames(Pan_TME[,-c(1:16)]), value = T))][x+1]),
               coef = cor.test(Pan_TME[,c(hallmark, grep(paste0(hallmark, "_"), colnames(Pan_TME[,-c(1:16)]), value = T))][,x+1],
                               Pan_TME[,c(hallmark, grep(paste0(hallmark, "_"), colnames(Pan_TME[,-c(1:16)]), value = T))][,1])$estimate,
               pval = cor.test(Pan_TME[,c(hallmark, grep(paste0(hallmark, "_"), colnames(Pan_TME[,-c(1:16)]), value = T))][,x+1],
                               Pan_TME[,c(hallmark, grep(paste0(hallmark, "_"), colnames(Pan_TME[,-c(1:16)]), value = T))][,1])$p.value)
  })
  #include results to the list
  corr.TME[[hallmark]] <- data.frame(t(tmp))
}

#adjust pvalues by using false discovery rate method
for (hallmark in names(corr.TME)) {
  corr.TME[[hallmark]]$adjP <- p.adjust(as.numeric(corr.TME[[hallmark]]$pval), method = "fdr", n= nrow(corr.TME[[hallmark]]))
  corr.TME[[hallmark]]$hallmark <- hallmark
}
#bind all data frames of the list 
corr.TME <- do.call(rbind, corr.TME)
#set columns as numeric
corr.TME$coef <- as.numeric(corr.TME$coef)
corr.TME$pval <- as.numeric(corr.TME$pval)
#set column as string/character
corr.TME$pathway <- as.character(corr.TME$pathway)

## Cancer
#create a list of empty data frames for each hallmark
corr.cancer <- list(H1 = data.frame(), H2= data.frame(), H3 = data.frame(),
                    H4 = data.frame(), H5 = data.frame(), H6 = data.frame(),
                    H7 = data.frame(), H8 = data.frame(), H9 = data.frame(),
                    H10 = data.frame(), H11 = data.frame(), H12 = data.frame(), H13 = data.frame())

#for each hallmark add the correlation (coefficient and significance) of hallmark activity and pathway score 
for (hallmark in paste0("H", 1:13)) {
  tmp <- sapply(1:ncol(Pan_cancer[,c(hallmark, grep(paste0(hallmark, "_"), colnames(Pan_cancer[,-c(1:16)]), value = T))][,-1]), function(x) {
    data.frame(pathway = colnames(Pan_cancer[,c(hallmark, grep(paste0(hallmark, "_"), colnames(Pan_cancer[,-c(1:16)]), value = T))][x+1]),
               coef = cor.test(Pan_cancer[,c(hallmark, grep(paste0(hallmark, "_"), colnames(Pan_cancer[,-c(1:16)]), value = T))][,x+1],
                               Pan_cancer[,c(hallmark, grep(paste0(hallmark, "_"), colnames(Pan_cancer[,-c(1:16)]), value = T))][,1])$estimate,
               pval = cor.test(Pan_cancer[,c(hallmark, grep(paste0(hallmark, "_"), colnames(Pan_cancer[,-c(1:16)]), value = T))][,x+1],
                               Pan_cancer[,c(hallmark, grep(paste0(hallmark, "_"), colnames(Pan_cancer[,-c(1:16)]), value = T))][,1])$p.value)
  })
  #include results to the list
  corr.cancer[[hallmark]] <- data.frame(t(tmp))
}

#adjust pvalues by using false discovery rate method
for (hallmark in names(corr.cancer)) {
  corr.cancer[[hallmark]]$adjP <- p.adjust(as.numeric(corr.cancer[[hallmark]]$pval), method = "fdr", n= nrow(corr.cancer[[hallmark]]))
  corr.cancer[[hallmark]]$hallmark <- hallmark
}
#bind all data frames of the list
corr.cancer <- do.call(rbind, corr.cancer)
#set columns as numeric
corr.cancer$coef <- as.numeric(corr.cancer$coef)
corr.cancer$pval <- as.numeric(corr.cancer$pval)
#set column as string/character
corr.cancer$pathway <- as.character(corr.cancer$pathway)

##### Plotting ####
library(ggpubr)
library(colorspace)
color_codes <-  list("Sustaining Proliferative Signaling"= "#15CE59",
                     "Evading Growth Suppressors" = "#701717",
                     "Avoiding Immune Destruction" = "#CB3BBD",
                     "Enabling Replicative Immortality" = "#4969D4",
                     "Tumour-Promoting Inflammation" = "#E6880D",
                     "Activating Invasion and Metastasis" = "#000000",
                     "Inducing Angiogenesis" = "#EE0F16",
                     "Genome Instability and Mutation" = "#132892",
                     "Resisting Cell Death" = "#8E909B",
                     "Deregulating Cellular Energetics" = "#71189E",
                     "Senescent cells" = "#05F3EB",
                     "Nonmutational Epigenetic reprogramming" = "#890269",
                     "Unlocking Phenotypic Plasticity" = "#95641A")

hallmark_names <- list("H1" = "Sustaining Proliferative Signaling",
                       "H2" = "Evading Growth Suppressors",
                       "H3" = "Avoiding Immune Destruction",
                       "H4" = "Enabling Replicative Immortality",
                       "H5" = "Tumour-Promoting Inflammation",
                       "H6" = "Activating Invasion and Metastasis",
                       "H7" = "Inducing Angiogenesis",
                       "H8" = "Genome Instability and Mutation",
                       "H9" = "Resisting Cell Death",
                       "H10"  = "Deregulating Cellular Energetics",
                       "H11"  = "Senescent cells",
                       "H12"  = "Nonmutational Epigenetic reprogramming",
                       "H13"  = "Unlocking Phenotypic Plasticity")

## TME
for (hallmark in paste0("H", 1:13)) {
  corr.TME$hallmark[corr.TME$hallmark == hallmark] <- hallmark_names[[hallmark]]
}

corr.TME$color <- corr.TME$hallmark

for (hallmark in names(color_codes)) {
  corr.TME$color[corr.TME$color == hallmark] <- color_codes[[hallmark]]
}

pathway_names <- sapply(corr.TME$pathway,function(x) {
  strsplit(x, split = "_")[[1]][2] })

corr.TME$pathway <- pathway_names

corr.TME$compartment <- "TME"

corr.TME.hallmarks <- filter(corr.TME, hallmark %in% c("Avoiding Immune Destruction", "Sustaining Proliferative Signaling",
                                                       "Resisting Cell Death",  "Inducing Angiogenesis", "Activating Invasion and Metastasis", "Tumour-Promoting Inflammation",  "Unlocking Phenotypic Plasticity"))

#corr.TME.hallmarks <- filter(corr.TME.hallmarks, coef > 0.5)

corr.TME.hallmarks <- arrange(corr.TME.hallmarks,hallmark, coef)
corr.TME.hallmarks$compartment <- "TME"


palette <- corr.TME.hallmarks[,c("hallmark", "color")]
palette <- data.frame(unique(palette), row.names = 1)

ggbarplot(corr.TME.hallmarks, x = "pathway", y = "coef",
          color = "hallmark", fill = "hallmark",
          palette = palette$color,
          label = F, lab.pos = "in", lab.col = "white",
          ggtheme = theme_pubclean(),
          #x.text.angle = 70,
          sort.by.groups = T,
          rotate = T
) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  font("y.text", size = 11) +
  font("title", color = darken("#00AFBB", amount = 0.3)) +
  font("caption", face = "italic") +
  font("legend.text", size = 13)


## Cancer

for (hallmark in paste0("H", 1:13)) {
  corr.cancer$hallmark[corr.cancer$hallmark == hallmark] <- hallmark_names[[hallmark]]
}

corr.cancer$color <- corr.cancer$hallmark

for (hallmark in names(color_codes)) {
  corr.cancer$color[corr.cancer$color == hallmark] <- color_codes[[hallmark]]
}

pathway_names <- sapply(corr.cancer$pathway,function(x) {
  strsplit(x, split = "_")[[1]][2] })

corr.cancer$pathway <- pathway_names

corr.cancer$compartment <- "Cancer"

corr.cancer.hallmarks <- filter(corr.cancer, !hallmark %in% c("Avoiding Immune Destruction", "Sustaining Proliferative Signaling",
                                                              "Resisting Cell Death",  "Inducing Angiogenesis", "Activating Invasion and Metastasis", "Tumour-Promoting Inflammation",  "Unlocking Phenotypic Plasticity"))

#corr.cancer.hallmarks <- filter(corr.cancer.hallmarks, coef > 0.5)

corr.cancer.hallmarks <- arrange(corr.cancer.hallmarks,hallmark, coef)
corr.cancer.hallmarks$compartment <- "Cancer"

palette <- corr.cancer.hallmarks[,c("hallmark", "color")]
palette <- data.frame(unique(palette), row.names = 1)

ggbarplot(corr.cancer.hallmarks, x = "pathway", y = "coef",
          color = "hallmark", fill = "hallmark",
          palette = palette$color,
          label = F, lab.pos = "in", lab.col = "white",
          ggtheme = theme_pubclean(),
          #x.text.angle = 70,
          sort.by.groups = T,
          rotate = T
) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  font("y.text", size = 11) +
  font("title", color = darken("#00AFBB", amount = 0.3)) +
  font("caption", face = "italic") +
  font("legend.text", size = 13)

## combine
corr.all <- rbind(corr.cancer.hallmarks, corr.TME.hallmarks)
corr.all <- filter(corr.all, coef > 0.5)

corr.all$hallmark <- factor(corr.all$hallmark, levels = c("Nonmutational Epigenetic reprogramming", "Enabling Replicative Immortality",
                                                          "Genome Instability and Mutation", "Evading Growth Suppressors", "Senescent cells",
                                                          "Deregulating Cellular Energetics", "Avoiding Immune Destruction", "Sustaining Proliferative Signaling",
                                                          "Resisting Cell Death",  "Inducing Angiogenesis", "Activating Invasion and Metastasis", "Tumour-Promoting Inflammation",  "Unlocking Phenotypic Plasticity"))
corr.all <- arrange(corr.all, hallmark, coef)

palette <- corr.all[,c("hallmark", "color")]
palette <- data.frame(unique(palette), row.names = 1)
#palette$color <- factor(palette$color, levels = palette$color)
#palette <- arrange(palette)

ggbarplot(corr.all, x = "pathway", y = "coef",
          color = "hallmark", fill = "hallmark",
          palette = palette$color,
          label = F, lab.pos = "in", lab.col = "white",
          ggtheme = theme_pubclean(),
          #y.text.angle = -45,
          facet.by = "compartment",
          sort.by.groups = T,
          rotate = T
) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  font("y.text", size = 11) +
  font("title", color = darken("#00AFBB", amount = 0.3)) +
  font("caption", face = "italic") +
  font("legend.text", size = 8)

### Grab top 3 pathways
corr.all <- rbind(corr.cancer.hallmarks, corr.TME.hallmarks)
corr.all <- filter(corr.all, coef > 0)

corr.all$hallmark <- factor(corr.all$hallmark, levels = c("Nonmutational Epigenetic reprogramming", "Enabling Replicative Immortality",
                                                          "Genome Instability and Mutation", "Evading Growth Suppressors", "Senescent cells",
                                                          "Deregulating Cellular Energetics", "Avoiding Immune Destruction", "Sustaining Proliferative Signaling",
                                                          "Resisting Cell Death",  "Inducing Angiogenesis", "Activating Invasion and Metastasis", "Tumour-Promoting Inflammation",  "Unlocking Phenotypic Plasticity"))
corr.all <- arrange(corr.all, hallmark, coef)

corr.all.top3 <- corr.all %>% dplyr::group_by(hallmark, compartment) %>% slice_max(n =3, order_by = coef, with_ties = F)

corr.all.top3 <- arrange(corr.all.top3, hallmark, coef)

palette <- corr.all.top3[,c("hallmark", "color")]
palette <- data.frame(unique(palette), row.names = 1)


ggbarplot(corr.all.top3, x = "pathway", y = "coef",
          color = "hallmark", fill = "hallmark",
          palette = palette$color,
          label = F, lab.pos = "in", lab.col = "white",
          ggtheme = theme_pubclean(),
          #y.text.angle = -45,
          facet.by = "compartment",
          sort.by.groups = T,
          rotate = T
) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  font("y.text", size = 15) +
  font("title", color = darken("#00AFBB", amount = 0.3)) +
  font("caption", face = "italic") +
  font("legend.text", size = 6)