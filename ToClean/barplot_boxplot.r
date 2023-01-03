annotation_method <- list("Breast" = "OCT",
                          "Ductal" = "OCT",
                          "IC" = "FFPE",
                          "Ovarian" = "OCT",
                          "Acinar" = "FFPE",
                          "Colorectal" = "OCT",
                          "Glioblastoma" = "OCT",
                          "OV4A" = "OCT",
                          #"Cervical" = "Cervix (n=1)",
                          "Intestine" = "OCT",
                          "OVFFPE" = "FFPE",
                          "DuctalFFPE" = "FFPE",
                          "ICC1L" = "OCT",
                          "cHC1T" = "OCT",
                          "HCC1T" ="OCT",
                          "BreastA" ="OCT",
                          "TNBCA" ="OCT",
                          "CUP295" ="FFPE",
                          "CUP297" = "FFPE",
                          "CUP298" = "FFPE",
                          "CUP303" = "FFPE",
                          "CUP306" = "FFPE",
                          "CUP309" = "FFPE",
                          "CUP319" = "FFPE",
                          "CUP312" = "FFPE",
                          "CUP327" = "FFPE",
                          "CUP114B" = "FFPE",
                          "OVD1" = "FFPE",
                          "HCC2T" = "OCT",
                          "CRC1" = "OCT",
                          "HCC5D" = "OCT",
                          "CRC2" = "OCT",
                          "DU1" = "FFPE",
                          "DU2" = "FFPE",
                          "DU3" = "FFPE",
                          "DU4" = "FFPE",
                          "DU5" = "FFPE",
                          "DU6" = "FFPE",
                          "DU7" = "FFPE",
                          "DU8" = "FFPE",
                          "DU9" = "FFPE",
                          "DU10" = "FFPE",
                          "DU11" = "FFPE",
                          "DU12" = "FFPE",
                          "DU13" = "FFPE",
                          "DU14" = "FFPE",
                          "DU15" = "FFPE",
                          "DU16" = "FFPE",
                          "C2" = "FFPE",
                          "C7" = "FFPE",
                          "C20" = "FFPE",
                          "C21" = "FFPE",
                          "C34" = "FFPE",
                          "C39" = "FFPE",
                          "C45" = "FFPE",
                          "C51" = "FFPE",
                          "frozen_A1" = "OCT",
                          "frozen_A3" = "OCT",
                          "frozen_C5" = "OCT",
                          "frozen_C23" = "OCT",
                          "P259_H2A2" = "OCT",
                          "P264" = "OCT",
                          "P270" = "OCT",
                          "P288" = "OCT",
                          "P306" = "OCT",
                          "UKF242T" = "OCT",
                          "UKF260T" = "OCT",
                          "UKF269T" = "OCT",
                          "UKF275T" = "OCT",
                          "PC1" = "OCT",
                          "PC2"= "OCT",
                          "P3" = "FFPE",
                          "P4" = "FFPE",
                          "P6" = "FFPE",
                          "P1" = "FFPE",
                          "P5" = "FFPE",
                          "P8" = "FFPE",
                          "M1" = "FFPE",
                          "M2" = "FFPE",
                          "M3" = "FFPE",
                          "M4" = "FFPE",
                          "Co1" = "FFPE",
                          "Co2" = "FFPE",
                          "Co3" = "FFPE",
                          "Co4" = "FFPE"
                          
)


### plot1: whole level

# library
library(tidyverse)
library(dplyr)
library(Seurat)

lm_mat.basal$diff <- lm_mat.1$Rsquared - lm_mat.basal$Rsquared
lm_mat.basal$RsquaredCancer <- lm_mat.basal.Cancer$Rsquared
lm_mat.basal$diffCancer <- lm_mat.1.Cancer$Rsquared - lm_mat.basal.Cancer$Rsquared
lm_mat.basal$RsquaredTME <- lm_mat.basal.TME$Rsquared
lm_mat.basal$diffTME <- lm_mat.1.TME$Rsquared - lm_mat.basal.TME$Rsquared

lm_mat.basal$method <- NA
for (sample in unique(lm_mat.basal$Sample)) {
  lm_mat.basal[lm_mat.basal$Sample==sample, "method"] <- annotation_method[[sample]]
}

#out negative values as 0 
lm_mat.basal[lm_mat.basal$diff < 0, "diff"] <- 0

lm_mat.basal[lm_mat.basal$diffCancer < 0, "diffCancer"] <- 0

lm_mat.basal[lm_mat.basal$diffTME < 0, "diffTME"] <- 0

lm_mat.basal[lm_mat.basal$diffBuffer < 0, "diffBuffer"] <- 0

df1 <- data.frame()

df1 <- rbind(df1, c(mean(lm_mat.basal$diff), "increase", "Whole"),
             c(mean(lm_mat.basal$Rsquared), "basal", "Whole"),
             c(mean(lm_mat.basal$diffCancer), "increase", "Cancer"),
             c(mean(lm_mat.basal$RsquaredCancer), "basal", "Cancer"),
             c(mean(lm_mat.basal$diffTME), "increase", "TME"),
             c(mean(lm_mat.basal$RsquaredTME), "basal", "TME"),
             c(mean(lm_mat.basal$diffBuffer), "increase", "Buffer"),
             c(mean(lm_mat.basal$RsquaredBuffer), "basal", "Buffer")
             )

colnames(df1) <- c("value", "model", "compartment")
df1$model <- factor(df1$model, levels = c("increase", "basal"))
df1$compartment <- factor(df1$compartment, levels = c("Whole", "Cancer", "Buffer", "TME"))
ggplot(df1, aes(x=compartment, y=as.numeric(value), fill=model)) + geom_bar(stat="identity")  +
  theme_classic() + labs(y="Mean R-squared", x="Tissue compartment")   + 
  theme(axis.ticks.x = element_blank()) + scale_fill_manual(values = c("darkred", "lightblue")) + 
  ylim(c(0, 1))

df1 <- rbind(df1, c(mean(lm_mat.basal$diff[lm_mat.basal$method=="OCT"]), "increase", "whole", "OCT"),
             c(mean(lm_mat.basal$Rsquared[lm_mat.basal$method=="OCT"]), "basal", "whole", "OCT"),
             c(mean(lm_mat.basal$diff[lm_mat.basal$method=="FFPE"]), "increase", "whole", "FFPE"),
             c(mean(lm_mat.basal$Rsquared[lm_mat.basal$method=="FFPE"]), "basal", "whole", "FFPE"),
             c(mean(lm_mat.basal$diffCancer[lm_mat.basal$method=="OCT"]), "increase", "Cancer", "OCT"), 
             c(mean(lm_mat.basal$RsquaredCancer[lm_mat.basal$method=="OCT"]), "basal", "Cancer", "OCT"),
             c(mean(lm_mat.basal$diffCancer[lm_mat.basal$method=="FFPE"]), "increase", "Cancer", "FFPE"), 
             c(mean(lm_mat.basal$RsquaredCancer[lm_mat.basal$method=="FFPE"]), "basal", "Cancer", "FFPE"),
             c(mean(lm_mat.basal$diffTME[lm_mat.basal$method=="OCT"]), "increase", "TME", "OCT"),
             c(mean(lm_mat.basal$RsquaredTME[lm_mat.basal$method=="OCT"]), "basal", "TME", "OCT"),
             c(mean(lm_mat.basal$diffTME[lm_mat.basal$method=="FFPE"]), "increase", "TME", "FFPE"),
             c(mean(lm_mat.basal$RsquaredTME[lm_mat.basal$method=="FFPE"]), "basal", "TME", "FFPE"),
             c(mean(lm_mat.basal$diffBuffer[lm_mat.basal$method=="OCT"]), "increase", "Buffer", "OCT"),
             c(mean(lm_mat.basal$RsquaredBuffer[lm_mat.basal$method=="OCT"]), "basal", "Buffer", "OCT"),
             c(mean(lm_mat.basal$diffBuffer[lm_mat.basal$method=="FFPE"]), "increase", "Buffer", "FFPE"),
             c(mean(lm_mat.basal$RsquaredBuffer[lm_mat.basal$method=="FFPE"]), "basal", "Buffer", "FFPE"))

colnames(df1) <- c("value", "model", "compartment", "method")
df1$id <- paste0(df1$compartment, "_", df1$method) 
df1$method[df1$model=="increase"] <- NA
df1$model <- factor(df1$model, levels = c("increase", "basal"))
df1$id <- factor(df1$id, levels = c("whole_OCT", "whole_FFPE" ,  "Cancer_OCT", "Cancer_FFPE", "TME_OCT", "TME_FFPE" ,  "Buffer_OCT", "Buffer_FFPE"))
ggplot(df1, aes(x=id, y=as.numeric(value), fill=model)) + geom_bar(stat="identity") + geom_text(aes(label=method, y=0.4)) +
  theme_classic() + labs(y="Mean R-squared", x="Tissue compartment")  + scale_x_discrete(labels=c("", "", "", "", "", "", "", "")) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(hjust=-1)) + scale_fill_manual(values = c("darkred", "lightblue"))
  


df2 <- final_data
df2 <- df2[!is.na(df2$y) & !df2$model %in% c("none_whole", "none_tme", "none_cancer"), ]
df2$facet <- "basal"
df2[df2$model %in% c("increase_whole", "increase_tme", "increase_cancer"), "facet"] <- "increase"
ggplot(df2, aes(x=model, y=y)) + geom_boxplot() + facet_wrap(~facet, scales = "free")


