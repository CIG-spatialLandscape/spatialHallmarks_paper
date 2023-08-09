##################################################
## Project: Cancer Hallmarks
## Script purpose: Establish control gene signatures based on two different perturbations
## Author: Sergi Cervilla* & Mustafa Sibai*
##################################################


library(stringr)
library(dplyr)
library(ggplot2)
# Import the list of pathways and genes associated with each Hallmark
paths_final <- read.delim("./paths_final.txt")
H.genes <- read.delim("./H.genes.txt")

####################### Perturbation 1: Randomly dropping genes per pathway ##########################

######### drop 25 % of genes randomly from each pathway and repeat 5 times #########
## R1
genes_R1 <- data.frame()

genes_R1 <- lapply(1:nrow(paths_final), function(x){
  rbind(genes_R1, data.frame(sample(unlist(str_split(paths_final$gene_list[x],", ")), size = round(length(unlist(str_split(paths_final$gene_list[x],", "))) *0.75))))
})

genes_R1 <- do.call(rbind, genes_R1)
colnames(genes_R1) <- "genes"

H.genes.R1 <- H.genes[H.genes$genes %in% genes_R1$genes,]


# Percentage of excluded genes from all of the list
(1- (length(unique(H.genes.R1$genes))/length(unique(H.genes$genes)))) * 100

# Percentage of excluded genes per Hallmark
(1- (table(H.genes.R1$H)/table(H.genes$H))) * 100

## R2
genes_R2 <- data.frame()

genes_R2 <- lapply(1:nrow(paths_final), function(x){
  rbind(genes_R2, data.frame(sample(unlist(str_split(paths_final$gene_list[x],", ")), size = round(length(unlist(str_split(paths_final$gene_list[x],", "))) *0.75))))
})

genes_R2 <- do.call(rbind, genes_R2)
colnames(genes_R2) <- "genes"

H.genes.R2 <- H.genes[H.genes$genes %in% genes_R2$genes,]


# Percentage of excluded genes from all of the list
(1- (length(unique(H.genes.R2$genes))/length(unique(H.genes$genes)))) * 100

# Percentage of excluded genes per Hallmark
(1- (table(H.genes.R2$H)/table(H.genes$H))) * 100

## R3
genes_R3 <- data.frame()

genes_R3 <- lapply(1:nrow(paths_final), function(x){
  rbind(genes_R3, data.frame(sample(unlist(str_split(paths_final$gene_list[x],", ")), size = round(length(unlist(str_split(paths_final$gene_list[x],", "))) *0.75))))
})

genes_R3 <- do.call(rbind, genes_R3)
colnames(genes_R3) <- "genes"

H.genes.R3 <- H.genes[H.genes$genes %in% genes_R3$genes,]


# Percentage of excluded genes from all of the list
(1- (length(unique(H.genes.R3$genes))/length(unique(H.genes$genes)))) * 100

# Percentage of excluded genes per Hallmark
(1- (table(H.genes.R3$H)/table(H.genes$H))) * 100


## R4
genes_R4 <- data.frame()

genes_R4 <- lapply(1:nrow(paths_final), function(x){
  rbind(genes_R4, data.frame(sample(unlist(str_split(paths_final$gene_list[x],", ")), size = round(length(unlist(str_split(paths_final$gene_list[x],", "))) *0.75))))
})

genes_R4 <- do.call(rbind, genes_R4)
colnames(genes_R4) <- "genes"

H.genes.R4 <- H.genes[H.genes$genes %in% genes_R4$genes,]


# Percentage of excluded genes from all of the list
(1- (length(unique(H.genes.R4$genes))/length(unique(H.genes$genes)))) * 100

# Percentage of excluded genes per Hallmark
(1- (table(H.genes.R4$H)/table(H.genes$H))) * 100

## R5
genes_R5 <- data.frame()

genes_R5 <- lapply(1:nrow(paths_final), function(x){
  rbind(genes_R5, data.frame(sample(unlist(str_split(paths_final$gene_list[x],", ")), size = round(length(unlist(str_split(paths_final$gene_list[x],", "))) *0.75))))
})

genes_R5 <- do.call(rbind, genes_R5)
colnames(genes_R5) <- "genes"

H.genes.R5 <- H.genes[H.genes$genes %in% genes_R5$genes,]


# Percentage of excluded genes from all of the list
(1- (length(unique(H.genes.R5$genes))/length(unique(H.genes$genes)))) * 100

# Percentage of excluded genes per Hallmark
(1- (table(H.genes.R5$H)/table(H.genes$H))) * 100



######### drop 50 % of genes randomly from each pathway and repeat 5 times #########
## R6
genes_R6 <- data.frame()

genes_R6 <- lapply(1:nrow(paths_final), function(x){
  rbind(genes_R6, data.frame(sample(unlist(str_split(paths_final$gene_list[x],", ")), size = round(length(unlist(str_split(paths_final$gene_list[x],", "))) *0.5))))
})

genes_R6 <- do.call(rbind, genes_R6)
colnames(genes_R6) <- "genes"

H.genes.R6 <- H.genes[H.genes$genes %in% genes_R6$genes,]


# Percentage of excluded genes from all of the list
(1- (length(unique(H.genes.R6$genes))/length(unique(H.genes$genes)))) * 100

# Percentage of excluded genes per Hallmark
(1- (table(H.genes.R6$H)/table(H.genes$H))) * 100


## R7
genes_R7 <- data.frame()

genes_R7 <- lapply(1:nrow(paths_final), function(x){
  rbind(genes_R7, data.frame(sample(unlist(str_split(paths_final$gene_list[x],", ")), size = round(length(unlist(str_split(paths_final$gene_list[x],", "))) *0.5))))
})

genes_R7 <- do.call(rbind, genes_R7)
colnames(genes_R7) <- "genes"

H.genes.R7 <- H.genes[H.genes$genes %in% genes_R7$genes,]


# Percentage of excluded genes from all of the list
(1- (length(unique(H.genes.R7$genes))/length(unique(H.genes$genes)))) * 100

# Percentage of excluded genes per Hallmark
(1- (table(H.genes.R7$H)/table(H.genes$H))) * 100

## R8
genes_R8 <- data.frame()

genes_R8 <- lapply(1:nrow(paths_final), function(x){
  rbind(genes_R8, data.frame(sample(unlist(str_split(paths_final$gene_list[x],", ")), size = round(length(unlist(str_split(paths_final$gene_list[x],", "))) *0.5))))
})

genes_R8 <- do.call(rbind, genes_R8)
colnames(genes_R8) <- "genes"

H.genes.R8 <- H.genes[H.genes$genes %in% genes_R8$genes,]


# Percentage of excluded genes from all of the list
(1- (length(unique(H.genes.R8$genes))/length(unique(H.genes$genes)))) * 100

# Percentage of excluded genes per Hallmark
(1- (table(H.genes.R8$H)/table(H.genes$H))) * 100

## R9
genes_R9 <- data.frame()

genes_R9 <- lapply(1:nrow(paths_final), function(x){
  rbind(genes_R9, data.frame(sample(unlist(str_split(paths_final$gene_list[x],", ")), size = round(length(unlist(str_split(paths_final$gene_list[x],", "))) *0.5))))
})

genes_R9 <- do.call(rbind, genes_R9)
colnames(genes_R9) <- "genes"

H.genes.R9 <- H.genes[H.genes$genes %in% genes_R9$genes,]


# Percentage of excluded genes from all of the list
(1- (length(unique(H.genes.R9$genes))/length(unique(H.genes$genes)))) * 100

# Percentage of excluded genes per Hallmark
(1- (table(H.genes.R9$H)/table(H.genes$H))) * 100

## R10
genes_R10 <- data.frame()

genes_R10 <- lapply(1:nrow(paths_final), function(x){
  rbind(genes_R10, data.frame(sample(unlist(str_split(paths_final$gene_list[x],", ")), size = round(length(unlist(str_split(paths_final$gene_list[x],", "))) *0.5))))
})

genes_R10 <- do.call(rbind, genes_R10)
colnames(genes_R10) <- "genes"

H.genes.R10 <- H.genes[H.genes$genes %in% genes_R10$genes,]


# Percentage of excluded genes from all of the list
(1- (length(unique(H.genes.R10$genes))/length(unique(H.genes$genes)))) * 100

# Percentage of excluded genes per Hallmark
(1- (table(H.genes.R10$H)/table(H.genes$H))) * 100


# R1
R1_pct <- data.frame((1- (table(H.genes.R1$H)/table(H.genes.gpt$H))) * 100)
colnames(R1_pct) <- c("Hallmark", "excl.pct")
R1_pct$exp <- "R1_0.25pct"
# R2
R2_pct <- data.frame((1- (table(H.genes.R2$H)/table(H.genes.gpt$H))) * 100)
colnames(R2_pct) <- c("Hallmark", "excl.pct")
R2_pct$exp <- "R2_0.25pct"
# R3
R3_pct <- data.frame((1- (table(H.genes.R3$H)/table(H.genes.gpt$H))) * 100)
colnames(R3_pct) <- c("Hallmark", "excl.pct")
R3_pct$exp <- "R3_0.25pct"
# R4
R4_pct <- data.frame((1- (table(H.genes.R4$H)/table(H.genes.gpt$H))) * 100)
colnames(R4_pct) <- c("Hallmark", "excl.pct")
R4_pct$exp <- "R4_0.25pct"
# R5
R5_pct <- data.frame((1- (table(H.genes.R5$H)/table(H.genes.gpt$H))) * 100)
colnames(R5_pct) <- c("Hallmark", "excl.pct")
R5_pct$exp <- "R5_0.25pct"
# R6
R6_pct <- data.frame((1- (table(H.genes.R6$H)/table(H.genes.gpt$H))) * 100)
colnames(R6_pct) <- c("Hallmark", "excl.pct")
R6_pct$exp <- "R6_0.50pct"
# R7
R7_pct <- data.frame((1- (table(H.genes.R7$H)/table(H.genes.gpt$H))) * 100)
colnames(R7_pct) <- c("Hallmark", "excl.pct")
R7_pct$exp <- "R7_0.50pct"
# R8
R8_pct <- data.frame((1- (table(H.genes.R8$H)/table(H.genes.gpt$H))) * 100)
colnames(R8_pct) <- c("Hallmark", "excl.pct")
R8_pct$exp <- "R8_0.50pct"
# R9
R9_pct <- data.frame((1- (table(H.genes.R9$H)/table(H.genes.gpt$H))) * 100)
colnames(R9_pct) <- c("Hallmark", "excl.pct")
R9_pct$exp <- "R9_0.50pct"
# R10
R10_pct <- data.frame((1- (table(H.genes.R10$H)/table(H.genes.gpt$H))) * 100)
colnames(R10_pct) <- c("Hallmark", "excl.pct")
R10_pct$exp <- "R10_0.50pct"

R_pct <- rbind(R1_pct, R2_pct, R3_pct, R4_pct, R5_pct, R6_pct, R7_pct, R8_pct, R9_pct, R10_pct)

R_pct$exp <- factor(R_pct$exp, levels= unique(R_pct$exp))

color_codes <-  c("H1"= "#15CE59", 
                  "H2" = "#701717",
                  "H3" = "#CB3BBD",
                  "H4" = "#4969D4",
                  "H5" = "#E6880D",
                  "H6" = "#000000",
                  "H7" = "#EE0F16",
                  "H8" = "#132892",
                  "H9" = "#8E909B",
                  "H10" = "#71189E",
                  "H11" = "#05F3EB",
                  "H12" = "#890269",
                  "H13" = "#95641A")
ggplot (R_pct,aes (x = Hallmark, fill = Hallmark, y = excl.pct)) + geom_bar(stat = "identity", position = "stack") + facet_wrap(~exp, nrow = 2) + scale_fill_manual(values = color_codes) +theme_classic() + ylim(c(0,100))



####################### Perturbation 2: Adding potentially important genes ##########################
##### get the overlap percentages of genes from each pathway of the filtered database with each hallmark gene set

# Paths_sum, H1.genes - H13.genes are from GeneCollection.R script

# H1
H1_intersect.pct <- data.frame(PATHWAY_NAMES = character(0), gene_list = character(0), Percentage = numeric(0), stringsAsFactors = FALSE)

for (i in seq_along(Paths_sum$gene_list)) {
  gene_list <- unlist(strsplit(Paths_sum$gene_list[i], ", "))
  intersect_genes <- intersect(H1.genes$genes, gene_list)
  overlap_count <- length(intersect_genes)
  percentage <- (overlap_count / length(gene_list)) * 100
  element_name <- Paths_sum$PATHWAY_NAMES[i]  # Use the pathway name from the 'pathway_name' column
  H1_intersect.pct <- rbind(H1_intersect.pct, data.frame(PATHWAY_NAMES = element_name, gene_list = paste(gene_list, collapse = ", "), Percentage = percentage, stringsAsFactors = FALSE))
}

H1_intersect.pct <- filter(H1_intersect.pct, Percentage >= 60 & Percentage <= 90)


H1_genes_unused <- data.frame(unlist(str_split(H1_intersect.pct$gene_list,", ")))
H1_genes_unused <- data.frame(H1_genes_unused[duplicated(H1_genes_unused) == FALSE,])
colnames(H1_genes_unused) <- "genes"
table(H1_genes_unused$genes %in% H1.genes$genes)
H1_genes_unused_unique <- data.frame(genes = H1_genes_unused[!H1_genes_unused$genes %in% H1.genes$genes,], H = "H1")



# H2
H2_intersect.pct <- data.frame(PATHWAY_NAMES = character(0), gene_list = character(0), Percentage = numeric(0), stringsAsFactors = FALSE)

for (i in seq_along(Paths_sum$gene_list)) {
  gene_list <- unlist(strsplit(Paths_sum$gene_list[i], ", "))
  intersect_genes <- intersect(H2.genes$genes, gene_list)
  overlap_count <- length(intersect_genes)
  percentage <- (overlap_count / length(gene_list)) * 100
  element_name <- Paths_sum$PATHWAY_NAMES[i]  # Use the pathway name from the 'pathway_name' column
  H2_intersect.pct <- rbind(H2_intersect.pct, data.frame(PATHWAY_NAMES = element_name, gene_list = paste(gene_list, collapse = ", "), Percentage = percentage, stringsAsFactors = FALSE))
}


H2_intersect.pct <- filter(H2_intersect.pct, Percentage >= 60 & Percentage <= 90)


H2_genes_unused <- data.frame(unlist(str_split(H2_intersect.pct$gene_list,", ")))
H2_genes_unused <- data.frame(H2_genes_unused[duplicated(H2_genes_unused) == FALSE,])
colnames(H2_genes_unused) <- "genes"
table(H2_genes_unused$genes %in% H2.genes$genes)
H2_genes_unused_unique <- data.frame(genes = H2_genes_unused[!H2_genes_unused$genes %in% H2.genes$genes,], H = "H2")


# H3
H3_intersect.pct <- data.frame(PATHWAY_NAMES = character(0), gene_list = character(0), Percentage = numeric(0), stringsAsFactors = FALSE)

for (i in seq_along(Paths_sum$gene_list)) {
  gene_list <- unlist(strsplit(Paths_sum$gene_list[i], ", "))
  intersect_genes <- intersect(H3.genes$genes, gene_list)
  overlap_count <- length(intersect_genes)
  percentage <- (overlap_count / length(gene_list)) * 100
  element_name <- Paths_sum$PATHWAY_NAMES[i]  # Use the pathway name from the 'pathway_name' column
  H3_intersect.pct <- rbind(H3_intersect.pct, data.frame(PATHWAY_NAMES = element_name, gene_list = paste(gene_list, collapse = ", "), Percentage = percentage, stringsAsFactors = FALSE))
}


H3_intersect.pct <- filter(H3_intersect.pct, Percentage >= 60 & Percentage <= 90)


H3_genes_unused <- data.frame(unlist(str_split(H3_intersect.pct$gene_list,", ")))
H3_genes_unused <- data.frame(H3_genes_unused[duplicated(H3_genes_unused) == FALSE,])
colnames(H3_genes_unused) <- "genes"
table(H3_genes_unused$genes %in% H3.genes$genes)
H3_genes_unused_unique <- data.frame(genes = H3_genes_unused[!H3_genes_unused$genes %in% H3.genes$genes,], H = "H3")

# H4
H4_intersect.pct <- data.frame(PATHWAY_NAMES = character(0), gene_list = character(0), Percentage = numeric(0), stringsAsFactors = FALSE)

for (i in seq_along(Paths_sum$gene_list)) {
  gene_list <- unlist(strsplit(Paths_sum$gene_list[i], ", "))
  intersect_genes <- intersect(H4.genes$genes, gene_list)
  overlap_count <- length(intersect_genes)
  percentage <- (overlap_count / length(gene_list)) * 100
  element_name <- Paths_sum$PATHWAY_NAMES[i]  # Use the pathway name from the 'pathway_name' column
  H4_intersect.pct <- rbind(H4_intersect.pct, data.frame(PATHWAY_NAMES = element_name, gene_list = paste(gene_list, collapse = ", "), Percentage = percentage, stringsAsFactors = FALSE))
}


H4_intersect.pct <- filter(H4_intersect.pct, Percentage >= 60 & Percentage <= 90)


H4_genes_unused <- data.frame(unlist(str_split(H4_intersect.pct$gene_list,", ")))
H4_genes_unused <- data.frame(H4_genes_unused[duplicated(H4_genes_unused) == FALSE,])
colnames(H4_genes_unused) <- "genes"
table(H4_genes_unused$genes %in% H4.genes$genes)
H4_genes_unused_unique <- data.frame(genes = H4_genes_unused[!H4_genes_unused$genes %in% H4.genes$genes,], H = "H4")


# H5
H5_intersect.pct <- data.frame(PATHWAY_NAMES = character(0), gene_list = character(0), Percentage = numeric(0), stringsAsFactors = FALSE)

for (i in seq_along(Paths_sum$gene_list)) {
  gene_list <- unlist(strsplit(Paths_sum$gene_list[i], ", "))
  intersect_genes <- intersect(H5.genes$genes, gene_list)
  overlap_count <- length(intersect_genes)
  percentage <- (overlap_count / length(gene_list)) * 100
  element_name <- Paths_sum$PATHWAY_NAMES[i]  # Use the pathway name from the 'pathway_name' column
  H5_intersect.pct <- rbind(H5_intersect.pct, data.frame(PATHWAY_NAMES = element_name, gene_list = paste(gene_list, collapse = ", "), Percentage = percentage, stringsAsFactors = FALSE))
}


H5_intersect.pct <- filter(H5_intersect.pct, Percentage >= 60 & Percentage <= 90)


H5_genes_unused <- data.frame(unlist(str_split(H5_intersect.pct$gene_list,", ")))
H5_genes_unused <- data.frame(H5_genes_unused[duplicated(H5_genes_unused) == FALSE,])
colnames(H5_genes_unused) <- "genes"
table(H5_genes_unused$genes %in% H5.genes$genes)
H5_genes_unused_unique <- data.frame(genes = H5_genes_unused[!H5_genes_unused$genes %in% H5.genes$genes,], H = "H5")

# H6
H6_intersect.pct <- data.frame(PATHWAY_NAMES = character(0), gene_list = character(0), Percentage = numeric(0), stringsAsFactors = FALSE)

for (i in seq_along(Paths_sum$gene_list)) {
  gene_list <- unlist(strsplit(Paths_sum$gene_list[i], ", "))
  intersect_genes <- intersect(H6.genes$genes, gene_list)
  overlap_count <- length(intersect_genes)
  percentage <- (overlap_count / length(gene_list)) * 100
  element_name <- Paths_sum$PATHWAY_NAMES[i]  # Use the pathway name from the 'pathway_name' column
  H6_intersect.pct <- rbind(H6_intersect.pct, data.frame(PATHWAY_NAMES = element_name, gene_list = paste(gene_list, collapse = ", "), Percentage = percentage, stringsAsFactors = FALSE))
}


H6_intersect.pct <- filter(H6_intersect.pct, Percentage >= 60 & Percentage <= 90)


H6_genes_unused <- data.frame(unlist(str_split(H6_intersect.pct$gene_list,", ")))
H6_genes_unused <- data.frame(H6_genes_unused[duplicated(H6_genes_unused) == FALSE,])
colnames(H6_genes_unused) <- "genes"
table(H6_genes_unused$genes %in% H6.genes$genes)
H6_genes_unused_unique <- data.frame(genes = H6_genes_unused[!H6_genes_unused$genes %in% H6.genes$genes,], H = "H6")

# H7
H7_intersect.pct <- data.frame(PATHWAY_NAMES = character(0), gene_list = character(0), Percentage = numeric(0), stringsAsFactors = FALSE)

for (i in seq_along(Paths_sum$gene_list)) {
  gene_list <- unlist(strsplit(Paths_sum$gene_list[i], ", "))
  intersect_genes <- intersect(H7.genes$genes, gene_list)
  overlap_count <- length(intersect_genes)
  percentage <- (overlap_count / length(gene_list)) * 100
  element_name <- Paths_sum$PATHWAY_NAMES[i]  # Use the pathway name from the 'pathway_name' column
  H7_intersect.pct <- rbind(H7_intersect.pct, data.frame(PATHWAY_NAMES = element_name, gene_list = paste(gene_list, collapse = ", "), Percentage = percentage, stringsAsFactors = FALSE))
}


H7_intersect.pct <- filter(H7_intersect.pct, Percentage >= 60 & Percentage <= 90)


H7_genes_unused <- data.frame(unlist(str_split(H7_intersect.pct$gene_list,", ")))
H7_genes_unused <- data.frame(H7_genes_unused[duplicated(H7_genes_unused) == FALSE,])
colnames(H7_genes_unused) <- "genes"
table(H7_genes_unused$genes %in% H7.genes$genes)
H7_genes_unused_unique <- data.frame(genes = H7_genes_unused[!H7_genes_unused$genes %in% H7.genes$genes,], H = "H7")

# H8
H8_intersect.pct <- data.frame(PATHWAY_NAMES = character(0), gene_list = character(0), Percentage = numeric(0), stringsAsFactors = FALSE)

for (i in seq_along(Paths_sum$gene_list)) {
  gene_list <- unlist(strsplit(Paths_sum$gene_list[i], ", "))
  intersect_genes <- intersect(H8.genes$genes, gene_list)
  overlap_count <- length(intersect_genes)
  percentage <- (overlap_count / length(gene_list)) * 100
  element_name <- Paths_sum$PATHWAY_NAMES[i]  # Use the pathway name from the 'pathway_name' column
  H8_intersect.pct <- rbind(H8_intersect.pct, data.frame(PATHWAY_NAMES = element_name, gene_list = paste(gene_list, collapse = ", "), Percentage = percentage, stringsAsFactors = FALSE))
}


H8_intersect.pct <- filter(H8_intersect.pct, Percentage >= 60 & Percentage <= 90)


H8_genes_unused <- data.frame(unlist(str_split(H8_intersect.pct$gene_list,", ")))
H8_genes_unused <- data.frame(H8_genes_unused[duplicated(H8_genes_unused) == FALSE,])
colnames(H8_genes_unused) <- "genes"
table(H8_genes_unused$genes %in% H8.genes$genes)
H8_genes_unused_unique <- data.frame(genes = H8_genes_unused[!H8_genes_unused$genes %in% H8.genes$genes,], H = "H8")

# H9
H9_intersect.pct <- data.frame(PATHWAY_NAMES = character(0), gene_list = character(0), Percentage = numeric(0), stringsAsFactors = FALSE)

for (i in seq_along(Paths_sum$gene_list)) {
  gene_list <- unlist(strsplit(Paths_sum$gene_list[i], ", "))
  intersect_genes <- intersect(H9.genes$genes, gene_list)
  overlap_count <- length(intersect_genes)
  percentage <- (overlap_count / length(gene_list)) * 100
  element_name <- Paths_sum$PATHWAY_NAMES[i]  # Use the pathway name from the 'pathway_name' column
  H9_intersect.pct <- rbind(H9_intersect.pct, data.frame(PATHWAY_NAMES = element_name, gene_list = paste(gene_list, collapse = ", "), Percentage = percentage, stringsAsFactors = FALSE))
}


H9_intersect.pct <- filter(H9_intersect.pct, Percentage >= 60 & Percentage <= 90)


H9_genes_unused <- data.frame(unlist(str_split(H9_intersect.pct$gene_list,", ")))
H9_genes_unused <- data.frame(H9_genes_unused[duplicated(H9_genes_unused) == FALSE,])
colnames(H9_genes_unused) <- "genes"
table(H9_genes_unused$genes %in% H9.genes$genes)
H9_genes_unused_unique <- data.frame(genes = H9_genes_unused[!H9_genes_unused$genes %in% H9.genes$genes,], H = "H9")


# H10
H10_intersect.pct <- data.frame(PATHWAY_NAMES = character(0), gene_list = character(0), Percentage = numeric(0), stringsAsFactors = FALSE)

for (i in seq_along(Paths_sum$gene_list)) {
  gene_list <- unlist(strsplit(Paths_sum$gene_list[i], ", "))
  intersect_genes <- intersect(H10.genes$genes, gene_list)
  overlap_count <- length(intersect_genes)
  percentage <- (overlap_count / length(gene_list)) * 100
  element_name <- Paths_sum$PATHWAY_NAMES[i]  # Use the pathway name from the 'pathway_name' column
  H10_intersect.pct <- rbind(H10_intersect.pct, data.frame(PATHWAY_NAMES = element_name, gene_list = paste(gene_list, collapse = ", "), Percentage = percentage, stringsAsFactors = FALSE))
}


H10_intersect.pct <- filter(H10_intersect.pct, Percentage >= 60 & Percentage <= 90)


H10_genes_unused <- data.frame(unlist(str_split(H10_intersect.pct$gene_list,", ")))
H10_genes_unused <- data.frame(H10_genes_unused[duplicated(H10_genes_unused) == FALSE,])
colnames(H10_genes_unused) <- "genes"
table(H10_genes_unused$genes %in% H10.genes$genes)
H10_genes_unused_unique <- data.frame(genes = H10_genes_unused[!H10_genes_unused$genes %in% H10.genes$genes,], H = "H10")

# H11
H11_intersect.pct <- data.frame(PATHWAY_NAMES = character(0), gene_list = character(0), Percentage = numeric(0), stringsAsFactors = FALSE)

for (i in seq_along(Paths_sum$gene_list)) {
  gene_list <- unlist(strsplit(Paths_sum$gene_list[i], ", "))
  intersect_genes <- intersect(H11.genes$genes, gene_list)
  overlap_count <- length(intersect_genes)
  percentage <- (overlap_count / length(gene_list)) * 100
  element_name <- Paths_sum$PATHWAY_NAMES[i]  # Use the pathway name from the 'pathway_name' column
  H11_intersect.pct <- rbind(H11_intersect.pct, data.frame(PATHWAY_NAMES = element_name, gene_list = paste(gene_list, collapse = ", "), Percentage = percentage, stringsAsFactors = FALSE))
}


H11_intersect.pct <- filter(H11_intersect.pct, Percentage >= 60 & Percentage <= 90)


H11_genes_unused <- data.frame(unlist(str_split(H11_intersect.pct$gene_list,", ")))
H11_genes_unused <- data.frame(H11_genes_unused[duplicated(H11_genes_unused) == FALSE,])
colnames(H11_genes_unused) <- "genes"
table(H11_genes_unused$genes %in% H11.genes$genes)
H11_genes_unused_unique <- data.frame(genes = H11_genes_unused[!H11_genes_unused$genes %in% H11.genes$genes,], H = "H11")

# H12
H12_intersect.pct <- data.frame(PATHWAY_NAMES = character(0), gene_list = character(0), Percentage = numeric(0), stringsAsFactors = FALSE)

for (i in seq_along(Paths_sum$gene_list)) {
  gene_list <- unlist(strsplit(Paths_sum$gene_list[i], ", "))
  intersect_genes <- intersect(H12.genes$genes, gene_list)
  overlap_count <- length(intersect_genes)
  percentage <- (overlap_count / length(gene_list)) * 100
  element_name <- Paths_sum$PATHWAY_NAMES[i]  # Use the pathway name from the 'pathway_name' column
  H12_intersect.pct <- rbind(H12_intersect.pct, data.frame(PATHWAY_NAMES = element_name, gene_list = paste(gene_list, collapse = ", "), Percentage = percentage, stringsAsFactors = FALSE))
}


H12_intersect.pct <- filter(H12_intersect.pct, Percentage >= 60 & Percentage <= 90)


H12_genes_unused <- data.frame(unlist(str_split(H12_intersect.pct$gene_list,", ")))
H12_genes_unused <- data.frame(H12_genes_unused[duplicated(H12_genes_unused) == FALSE,])
colnames(H12_genes_unused) <- "genes"
table(H12_genes_unused$genes %in% H12.genes$genes)
H12_genes_unused_unique <- data.frame(genes = H12_genes_unused[!H12_genes_unused$genes %in% H12.genes$genes,], H = "H12")

# H13
H13_intersect.pct <- data.frame(PATHWAY_NAMES = character(0), gene_list = character(0), Percentage = numeric(0), stringsAsFactors = FALSE)

for (i in seq_along(Paths_sum$gene_list)) {
  gene_list <- unlist(strsplit(Paths_sum$gene_list[i], ", "))
  intersect_genes <- intersect(H13.genes$genes, gene_list)
  overlap_count <- length(intersect_genes)
  percentage <- (overlap_count / length(gene_list)) * 100
  element_name <- Paths_sum$PATHWAY_NAMES[i]  # Use the pathway name from the 'pathway_name' column
  H13_intersect.pct <- rbind(H13_intersect.pct, data.frame(PATHWAY_NAMES = element_name, gene_list = paste(gene_list, collapse = ", "), Percentage = percentage, stringsAsFactors = FALSE))
}


H13_intersect.pct <- filter(H13_intersect.pct, Percentage >= 60 & Percentage <= 90)


H13_genes_unused <- data.frame(unlist(str_split(H13_intersect.pct$gene_list,", ")))
H13_genes_unused <- data.frame(H13_genes_unused[duplicated(H13_genes_unused) == FALSE,])
colnames(H13_genes_unused) <- "genes"
table(H13_genes_unused$genes %in% H13.genes$genes)
H13_genes_unused_unique <- data.frame(genes = H13_genes_unused[!H13_genes_unused$genes %in% H13.genes$genes,], H = "H13")


#################### pick a number of random unique genes corresponding to 5%, 10%, 25%, and 50% percent of the total size of the tested Hallmark's gene set. Apply this 5 times for each percentage. 
#################### Then, add these additional genes to the original gene set of the tested Hallmark #########

########### 5 pct

# H1
H1_unique_5pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H1_genes_unused_unique$genes, size = round(length(H1.genes$genes) *0.05)), H = "H1")
})

for (df in 1:length(H1_unique_5pct)){
  H1_unique_5pct[[df]] <- rbind(H1_unique_5pct[[df]], H1.genes)
}
# H2
H2_unique_5pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H2_genes_unused_unique$genes, size = round(length(H2.genes$genes) *0.05)), H = "H2")
})

for (df in 1:length(H2_unique_5pct)){
  H2_unique_5pct[[df]] <- rbind(H2_unique_5pct[[df]], H2.genes)
}

# H3
H3_unique_5pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H3_genes_unused_unique$genes, size = round(length(H3.genes$genes) *0.05)), H = "H3")
})

for (df in 1:length(H3_unique_5pct)){
  H3_unique_5pct[[df]] <- rbind(H3_unique_5pct[[df]], H3.genes)
}

# H4
H4_unique_5pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H4_genes_unused_unique$genes, size = round(length(H4.genes$genes) *0.05)), H = "H4")
})

for (df in 1:length(H4_unique_5pct)){
  H4_unique_5pct[[df]] <- rbind(H4_unique_5pct[[df]], H4.genes)
}

# H5
H5_unique_5pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H5_genes_unused_unique$genes, size = round(length(H5.genes$genes) *0.05)), H = "H5")
})

for (df in 1:length(H5_unique_5pct)){
  H5_unique_5pct[[df]] <- rbind(H5_unique_5pct[[df]], H5.genes)
}

# H6
H6_unique_5pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H6_genes_unused_unique$genes, size = round(length(H6.genes$genes) *0.05)), H = "H6")
})

for (df in 1:length(H6_unique_5pct)){
  H6_unique_5pct[[df]] <- rbind(H6_unique_5pct[[df]], H6.genes)
}

# H7
H7_unique_5pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H7_genes_unused_unique$genes, size = round(length(H7.genes$genes) *0.05)), H = "H7")
})

for (df in 1:length(H7_unique_5pct)){
  H7_unique_5pct[[df]] <- rbind(H7_unique_5pct[[df]], H7.genes)
}

# H8
H8_unique_5pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H8_genes_unused_unique$genes, size = round(length(H8.genes$genes) *0.05)), H = "H8")
})

for (df in 1:length(H8_unique_5pct)){
  H8_unique_5pct[[df]] <- rbind(H8_unique_5pct[[df]], H8.genes)
}

# H9
H9_unique_5pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H9_genes_unused_unique$genes, size = round(length(H9.genes$genes) *0.05)), H = "H9")
})

for (df in 1:length(H9_unique_5pct)){
  H9_unique_5pct[[df]] <- rbind(H9_unique_5pct[[df]], H9.genes)
}

# H10
H10_unique_5pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H10_genes_unused_unique$genes, size = round(length(H10.genes$genes) *0.05)), H = "H10")
})

for (df in 1:length(H10_unique_5pct)){
  H10_unique_5pct[[df]] <- rbind(H10_unique_5pct[[df]], H10.genes)
}

# H11
H11_unique_5pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H11_genes_unused_unique$genes, size = round(length(H11.genes$genes) *0.05)), H = "H11")
})

for (df in 1:length(H11_unique_5pct)){
  H11_unique_5pct[[df]] <- rbind(H11_unique_5pct[[df]], H11.genes)
}

# H12
H12_unique_5pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H12_genes_unused_unique$genes, size = round(length(H12.genes$genes) *0.05)), H = "H12")
})

for (df in 1:length(H12_unique_5pct)){
  H12_unique_5pct[[df]] <- rbind(H12_unique_5pct[[df]], H12.genes)
}

# H13
H13_unique_5pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H13_genes_unused_unique$genes, size = round(length(H13.genes$genes) *0.05)), H = "H13")
})

for (df in 1:length(H13_unique_5pct)){
  H13_unique_5pct[[df]] <- rbind(H13_unique_5pct[[df]], H13.genes)
}


H_genes_gpt_5pct_R1 <- rbind(H1_unique_5pct[[1]], H2_unique_5pct[[1]], H3_unique_5pct[[1]], H4_unique_5pct[[1]],
                             H5_unique_5pct[[1]], H6_unique_5pct[[1]], H7_unique_5pct[[1]], H8_unique_5pct[[1]],
                             H9_unique_5pct[[1]], H10_unique_5pct[[1]], H11_unique_5pct[[1]], H12_unique_5pct[[1]],
                             H13_unique_5pct[[1]])
H_genes_gpt_5pct_R2 <- rbind(H1_unique_5pct[[2]], H2_unique_5pct[[2]], H3_unique_5pct[[2]], H4_unique_5pct[[2]],
                             H5_unique_5pct[[2]], H6_unique_5pct[[2]], H7_unique_5pct[[2]], H8_unique_5pct[[2]],
                             H9_unique_5pct[[2]], H10_unique_5pct[[2]], H11_unique_5pct[[2]], H12_unique_5pct[[2]],
                             H13_unique_5pct[[2]])
H_genes_gpt_5pct_R3 <- rbind(H1_unique_5pct[[3]], H2_unique_5pct[[3]], H3_unique_5pct[[3]], H4_unique_5pct[[3]],
                             H5_unique_5pct[[3]], H6_unique_5pct[[3]], H7_unique_5pct[[3]], H8_unique_5pct[[3]],
                             H9_unique_5pct[[3]], H10_unique_5pct[[3]], H11_unique_5pct[[3]], H12_unique_5pct[[3]],
                             H13_unique_5pct[[3]])
H_genes_gpt_5pct_R4 <- rbind(H1_unique_5pct[[4]], H2_unique_5pct[[4]], H3_unique_5pct[[4]], H4_unique_5pct[[4]],
                             H5_unique_5pct[[4]], H6_unique_5pct[[4]], H7_unique_5pct[[4]], H8_unique_5pct[[4]],
                             H9_unique_5pct[[4]], H10_unique_5pct[[4]], H11_unique_5pct[[4]], H12_unique_5pct[[4]],
                             H13_unique_5pct[[4]])
H_genes_gpt_5pct_R5 <- rbind(H1_unique_5pct[[5]], H2_unique_5pct[[5]], H3_unique_5pct[[5]], H4_unique_5pct[[5]],
                             H5_unique_5pct[[5]], H6_unique_5pct[[5]], H7_unique_5pct[[5]], H8_unique_5pct[[5]],
                             H9_unique_5pct[[5]], H10_unique_5pct[[5]], H11_unique_5pct[[5]], H12_unique_5pct[[5]],
                             H13_unique_5pct[[5]])


####### 10 pct

# H1
H1_unique_10pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H1_genes_unused_unique$genes, size = round(length(H1.genes$genes) *0.10)), H = "H1")
})

for (df in 1:length(H1_unique_10pct)){
  H1_unique_10pct[[df]] <- rbind(H1_unique_10pct[[df]], H1.genes)
}
# H2
H2_unique_10pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H2_genes_unused_unique$genes, size = round(length(H2.genes$genes) *0.10)), H = "H2")
})

for (df in 1:length(H2_unique_10pct)){
  H2_unique_10pct[[df]] <- rbind(H2_unique_10pct[[df]], H2.genes)
}

# H3
H3_unique_10pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H3_genes_unused_unique$genes, size = round(length(H3.genes$genes) *0.10)), H = "H3")
})

for (df in 1:length(H3_unique_10pct)){
  H3_unique_10pct[[df]] <- rbind(H3_unique_10pct[[df]], H3.genes)
}

# H4
H4_unique_10pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H4_genes_unused_unique$genes, size = round(length(H4.genes$genes) *0.10)), H = "H4")
})

for (df in 1:length(H4_unique_10pct)){
  H4_unique_10pct[[df]] <- rbind(H4_unique_10pct[[df]], H4.genes)
}

# H5
H5_unique_10pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H5_genes_unused_unique$genes, size = round(length(H5.genes$genes) *0.10)), H = "H5")
})

for (df in 1:length(H5_unique_10pct)){
  H5_unique_10pct[[df]] <- rbind(H5_unique_10pct[[df]], H5.genes)
}

# H6
H6_unique_10pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H6_genes_unused_unique$genes, size = round(length(H6.genes$genes) *0.10)), H = "H6")
})

for (df in 1:length(H6_unique_10pct)){
  H6_unique_10pct[[df]] <- rbind(H6_unique_10pct[[df]], H6.genes)
}

# H7
H7_unique_10pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H7_genes_unused_unique$genes, size = round(length(H7.genes$genes) *0.10)), H = "H7")
})

for (df in 1:length(H7_unique_10pct)){
  H7_unique_10pct[[df]] <- rbind(H7_unique_10pct[[df]], H7.genes)
}

# H8
H8_unique_10pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H8_genes_unused_unique$genes, size = round(length(H8.genes$genes) *0.10)), H = "H8")
})

for (df in 1:length(H8_unique_10pct)){
  H8_unique_10pct[[df]] <- rbind(H8_unique_10pct[[df]], H8.genes)
}

# H9
H9_unique_10pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H9_genes_unused_unique$genes, size = round(length(H9.genes$genes) *0.10)), H = "H9")
})

for (df in 1:length(H9_unique_10pct)){
  H9_unique_10pct[[df]] <- rbind(H9_unique_10pct[[df]], H9.genes)
}

# H10
H10_unique_10pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H10_genes_unused_unique$genes, size = round(length(H10.genes$genes) *0.10)), H = "H10")
})

for (df in 1:length(H10_unique_10pct)){
  H10_unique_10pct[[df]] <- rbind(H10_unique_10pct[[df]], H10.genes)
}

# H11
H11_unique_10pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H11_genes_unused_unique$genes, size = round(length(H11.genes$genes) *0.10)), H = "H11")
})

for (df in 1:length(H11_unique_10pct)){
  H11_unique_10pct[[df]] <- rbind(H11_unique_10pct[[df]], H11.genes)
}

# H12
H12_unique_10pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H12_genes_unused_unique$genes, size = round(length(H12.genes$genes) *0.10)), H = "H12")
})

for (df in 1:length(H12_unique_10pct)){
  H12_unique_10pct[[df]] <- rbind(H12_unique_10pct[[df]], H12.genes)
}

# H13
H13_unique_10pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H13_genes_unused_unique$genes, size = round(length(H13.genes$genes) *0.10)), H = "H13")
})

for (df in 1:length(H13_unique_10pct)){
  H13_unique_10pct[[df]] <- rbind(H13_unique_10pct[[df]], H13.genes)
}



H_genes_gpt_10pct_R1 <- rbind(H1_unique_10pct[[1]], H2_unique_10pct[[1]], H3_unique_10pct[[1]], H4_unique_10pct[[1]],
                              H5_unique_10pct[[1]], H6_unique_10pct[[1]], H7_unique_10pct[[1]], H8_unique_10pct[[1]],
                              H9_unique_10pct[[1]], H10_unique_10pct[[1]], H11_unique_10pct[[1]], H12_unique_10pct[[1]],
                              H13_unique_10pct[[1]])
H_genes_gpt_10pct_R2 <- rbind(H1_unique_10pct[[2]], H2_unique_10pct[[2]], H3_unique_10pct[[2]], H4_unique_10pct[[2]],
                              H5_unique_10pct[[2]], H6_unique_10pct[[2]], H7_unique_10pct[[2]], H8_unique_10pct[[2]],
                              H9_unique_10pct[[2]], H10_unique_10pct[[2]], H11_unique_10pct[[2]], H12_unique_10pct[[2]],
                              H13_unique_10pct[[2]])
H_genes_gpt_10pct_R3 <- rbind(H1_unique_10pct[[3]], H2_unique_10pct[[3]], H3_unique_10pct[[3]], H4_unique_10pct[[3]],
                              H5_unique_10pct[[3]], H6_unique_10pct[[3]], H7_unique_10pct[[3]], H8_unique_10pct[[3]],
                              H9_unique_10pct[[3]], H10_unique_10pct[[3]], H11_unique_10pct[[3]], H12_unique_10pct[[3]],
                              H13_unique_10pct[[3]])
H_genes_gpt_10pct_R4 <- rbind(H1_unique_10pct[[4]], H2_unique_10pct[[4]], H3_unique_10pct[[4]], H4_unique_10pct[[4]],
                              H5_unique_10pct[[4]], H6_unique_10pct[[4]], H7_unique_10pct[[4]], H8_unique_10pct[[4]],
                              H9_unique_10pct[[4]], H10_unique_10pct[[4]], H11_unique_10pct[[4]], H12_unique_10pct[[4]],
                              H13_unique_10pct[[4]])
H_genes_gpt_10pct_R5 <- rbind(H1_unique_10pct[[5]], H2_unique_10pct[[5]], H3_unique_10pct[[5]], H4_unique_10pct[[5]],
                              H5_unique_10pct[[5]], H6_unique_10pct[[5]], H7_unique_10pct[[5]], H8_unique_10pct[[5]],
                              H9_unique_10pct[[5]], H10_unique_10pct[[5]], H11_unique_10pct[[5]], H12_unique_10pct[[5]],
                              H13_unique_10pct[[5]])


########## 25 pct

# H1
H1_unique_25pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H1_genes_unused_unique$genes, size = round(length(H1.genes$genes) *0.25)), H = "H1")
})

for (df in 1:length(H1_unique_25pct)){
  H1_unique_25pct[[df]] <- rbind(H1_unique_25pct[[df]], H1.genes)
}
# H2
H2_unique_25pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H2_genes_unused_unique$genes, size = round(length(H2.genes$genes) *0.25)), H = "H2")
})

for (df in 1:length(H2_unique_25pct)){
  H2_unique_25pct[[df]] <- rbind(H2_unique_25pct[[df]], H2.genes)
}

# H3
H3_unique_25pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H3_genes_unused_unique$genes, size = round(length(H3.genes$genes) *0.25)), H = "H3")
})

for (df in 1:length(H3_unique_25pct)){
  H3_unique_25pct[[df]] <- rbind(H3_unique_25pct[[df]], H3.genes)
}

# H4
H4_unique_25pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H4_genes_unused_unique$genes, size = round(length(H4.genes$genes) *0.25)), H = "H4")
})

for (df in 1:length(H4_unique_25pct)){
  H4_unique_25pct[[df]] <- rbind(H4_unique_25pct[[df]], H4.genes)
}

# H5
H5_unique_25pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H5_genes_unused_unique$genes, size = round(length(H5.genes$genes) *0.25)), H = "H5")
})

for (df in 1:length(H5_unique_25pct)){
  H5_unique_25pct[[df]] <- rbind(H5_unique_25pct[[df]], H5.genes)
}

# H6
H6_unique_25pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H6_genes_unused_unique$genes, size = round(length(H6.genes$genes) *0.25)), H = "H6")
})

for (df in 1:length(H6_unique_25pct)){
  H6_unique_25pct[[df]] <- rbind(H6_unique_25pct[[df]], H6.genes)
}

# H7
H7_unique_25pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H7_genes_unused_unique$genes, size = round(length(H7.genes$genes) *0.25)), H = "H7")
})

for (df in 1:length(H7_unique_25pct)){
  H7_unique_25pct[[df]] <- rbind(H7_unique_25pct[[df]], H7.genes)
}

# H8
H8_unique_25pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H8_genes_unused_unique$genes, size = round(length(H8.genes$genes) *0.25)), H = "H8")
})

for (df in 1:length(H8_unique_25pct)){
  H8_unique_25pct[[df]] <- rbind(H8_unique_25pct[[df]], H8.genes)
}

# H9
H9_unique_25pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H9_genes_unused_unique$genes, size = round(length(H9.genes$genes) *0.25)), H = "H9")
})

for (df in 1:length(H9_unique_25pct)){
  H9_unique_25pct[[df]] <- rbind(H9_unique_25pct[[df]], H9.genes)
}

# H10
H10_unique_25pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H10_genes_unused_unique$genes, size = round(length(H10.genes$genes) *0.25)), H = "H10")
})

for (df in 1:length(H10_unique_25pct)){
  H10_unique_25pct[[df]] <- rbind(H10_unique_25pct[[df]], H10.genes)
}

# H11
H11_unique_25pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H11_genes_unused_unique$genes, size = round(length(H11.genes$genes) *0.25)), H = "H11")
})

for (df in 1:length(H11_unique_25pct)){
  H11_unique_25pct[[df]] <- rbind(H11_unique_25pct[[df]], H11.genes)
}

# H12
H12_unique_25pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H12_genes_unused_unique$genes, size = round(length(H12.genes$genes) *0.25)), H = "H12")
})

for (df in 1:length(H12_unique_25pct)){
  H12_unique_25pct[[df]] <- rbind(H12_unique_25pct[[df]], H12.genes)
}

# H13
H13_unique_25pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H13_genes_unused_unique$genes, size = round(length(H13.genes$genes) *0.25)), H = "H13")
})

for (df in 1:length(H13_unique_25pct)){
  H13_unique_25pct[[df]] <- rbind(H13_unique_25pct[[df]], H13.genes)
}


H_genes_gpt_25pct_R1 <- rbind(H1_unique_25pct[[1]], H2_unique_25pct[[1]], H3_unique_25pct[[1]], H4_unique_25pct[[1]],
                              H5_unique_25pct[[1]], H6_unique_25pct[[1]], H7_unique_25pct[[1]], H8_unique_25pct[[1]],
                              H9_unique_25pct[[1]], H10_unique_25pct[[1]], H11_unique_25pct[[1]], H12_unique_25pct[[1]],
                              H13_unique_25pct[[1]])
H_genes_gpt_25pct_R2 <- rbind(H1_unique_25pct[[2]], H2_unique_25pct[[2]], H3_unique_25pct[[2]], H4_unique_25pct[[2]],
                              H5_unique_25pct[[2]], H6_unique_25pct[[2]], H7_unique_25pct[[2]], H8_unique_25pct[[2]],
                              H9_unique_25pct[[2]], H10_unique_25pct[[2]], H11_unique_25pct[[2]], H12_unique_25pct[[2]],
                              H13_unique_25pct[[2]])
H_genes_gpt_25pct_R3 <- rbind(H1_unique_25pct[[3]], H2_unique_25pct[[3]], H3_unique_25pct[[3]], H4_unique_25pct[[3]],
                              H5_unique_25pct[[3]], H6_unique_25pct[[3]], H7_unique_25pct[[3]], H8_unique_25pct[[3]],
                              H9_unique_25pct[[3]], H10_unique_25pct[[3]], H11_unique_25pct[[3]], H12_unique_25pct[[3]],
                              H13_unique_25pct[[3]])
H_genes_gpt_25pct_R4 <- rbind(H1_unique_25pct[[4]], H2_unique_25pct[[4]], H3_unique_25pct[[4]], H4_unique_25pct[[4]],
                              H5_unique_25pct[[4]], H6_unique_25pct[[4]], H7_unique_25pct[[4]], H8_unique_25pct[[4]],
                              H9_unique_25pct[[4]], H10_unique_25pct[[4]], H11_unique_25pct[[4]], H12_unique_25pct[[4]],
                              H13_unique_25pct[[4]])
H_genes_gpt_25pct_R5 <- rbind(H1_unique_25pct[[5]], H2_unique_25pct[[5]], H3_unique_25pct[[5]], H4_unique_25pct[[5]],
                              H5_unique_25pct[[5]], H6_unique_25pct[[5]], H7_unique_25pct[[5]], H8_unique_25pct[[5]],
                              H9_unique_25pct[[5]], H10_unique_25pct[[5]], H11_unique_25pct[[5]], H12_unique_25pct[[5]],
                              H13_unique_25pct[[5]])



######### 50 pct

# H1
H1_unique_50pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H1_genes_unused_unique$genes, size = round(length(H1.genes$genes) *0.50)), H = "H1")
})

for (df in 1:length(H1_unique_50pct)){
  H1_unique_50pct[[df]] <- rbind(H1_unique_50pct[[df]], H1.genes)
}
# H2
H2_unique_50pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H2_genes_unused_unique$genes, size = round(length(H2.genes$genes) *0.50)), H = "H2")
})

for (df in 1:length(H2_unique_50pct)){
  H2_unique_50pct[[df]] <- rbind(H2_unique_50pct[[df]], H2.genes)
}

# H3
H3_unique_50pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H3_genes_unused_unique$genes, size = round(length(H3.genes$genes) *0.50)), H = "H3")
})

for (df in 1:length(H3_unique_50pct)){
  H3_unique_50pct[[df]] <- rbind(H3_unique_50pct[[df]], H3.genes)
}

# H4
H4_unique_50pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H4_genes_unused_unique$genes, size = round(length(H4.genes$genes) *0.50)), H = "H4")
})

for (df in 1:length(H4_unique_50pct)){
  H4_unique_50pct[[df]] <- rbind(H4_unique_50pct[[df]], H4.genes)
}

# H5
H5_unique_50pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H5_genes_unused_unique$genes, size = round(length(H5.genes$genes) *0.3)), H = "H5")
})

for (df in 1:length(H5_unique_50pct)){
  H5_unique_50pct[[df]] <- rbind(H5_unique_50pct[[df]], H5.genes)
}

# H6
H6_unique_50pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H6_genes_unused_unique$genes, size = round(length(H6.genes$genes) *0.50)), H = "H6")
})

for (df in 1:length(H6_unique_50pct)){
  H6_unique_50pct[[df]] <- rbind(H6_unique_50pct[[df]], H6.genes)
}

# H7
H7_unique_50pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H7_genes_unused_unique$genes, size = round(length(H7.genes$genes) *0.50)), H = "H7")
})

for (df in 1:length(H7_unique_50pct)){
  H7_unique_50pct[[df]] <- rbind(H7_unique_50pct[[df]], H7.genes)
}

# H8
H8_unique_50pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H8_genes_unused_unique$genes, size = round(length(H8.genes$genes) *0.50)), H = "H8")
})

for (df in 1:length(H8_unique_50pct)){
  H8_unique_50pct[[df]] <- rbind(H8_unique_50pct[[df]], H8.genes)
}

# H9
H9_unique_50pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H9_genes_unused_unique$genes, size = round(length(H9.genes$genes) *0.50)), H = "H9")
})

for (df in 1:length(H9_unique_50pct)){
  H9_unique_50pct[[df]] <- rbind(H9_unique_50pct[[df]], H9.genes)
}

# H10
H10_unique_50pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H10_genes_unused_unique$genes, size = round(length(H10.genes$genes) *0.50)), H = "H10")
})

for (df in 1:length(H10_unique_50pct)){
  H10_unique_50pct[[df]] <- rbind(H10_unique_50pct[[df]], H10.genes)
}

# H11
H11_unique_50pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H11_genes_unused_unique$genes, size = round(length(H11.genes$genes) *0.50)), H = "H11")
})

for (df in 1:length(H11_unique_50pct)){
  H11_unique_50pct[[df]] <- rbind(H11_unique_50pct[[df]], H11.genes)
}

# H12
H12_unique_50pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H12_genes_unused_unique$genes, size = round(length(H12.genes$genes) *0.50)), H = "H12")
})

for (df in 1:length(H12_unique_50pct)){
  H12_unique_50pct[[df]] <- rbind(H12_unique_50pct[[df]], H12.genes)
}

# H13
H13_unique_50pct <- lapply(1:5, function(x){
  data.frame(genes = sample(H13_genes_unused_unique$genes, size = round(length(H13.genes$genes) *0.50)), H = "H13")
})

for (df in 1:length(H13_unique_50pct)){
  H13_unique_50pct[[df]] <- rbind(H13_unique_50pct[[df]], H13.genes)
}


H_genes_gpt_50pct_R1 <- rbind(H1_unique_50pct[[1]], H2_unique_50pct[[1]], H3_unique_50pct[[1]], H4_unique_50pct[[1]],
                              H5_unique_50pct[[1]], H6_unique_50pct[[1]], H7_unique_50pct[[1]], H8_unique_50pct[[1]],
                              H9_unique_50pct[[1]], H10_unique_50pct[[1]], H11_unique_50pct[[1]], H12_unique_50pct[[1]],
                              H13_unique_50pct[[1]])
H_genes_gpt_50pct_R2 <- rbind(H1_unique_50pct[[2]], H2_unique_50pct[[2]], H3_unique_50pct[[2]], H4_unique_50pct[[2]],
                              H5_unique_50pct[[2]], H6_unique_50pct[[2]], H7_unique_50pct[[2]], H8_unique_50pct[[2]],
                              H9_unique_50pct[[2]], H10_unique_50pct[[2]], H11_unique_50pct[[2]], H12_unique_50pct[[2]],
                              H13_unique_50pct[[2]])
H_genes_gpt_50pct_R3 <- rbind(H1_unique_50pct[[3]], H2_unique_50pct[[3]], H3_unique_50pct[[3]], H4_unique_50pct[[3]],
                              H5_unique_50pct[[3]], H6_unique_50pct[[3]], H7_unique_50pct[[3]], H8_unique_50pct[[3]],
                              H9_unique_50pct[[3]], H10_unique_50pct[[3]], H11_unique_50pct[[3]], H12_unique_50pct[[3]],
                              H13_unique_50pct[[3]])
H_genes_gpt_50pct_R4 <- rbind(H1_unique_50pct[[4]], H2_unique_50pct[[4]], H3_unique_50pct[[4]], H4_unique_50pct[[4]],
                              H5_unique_50pct[[4]], H6_unique_50pct[[4]], H7_unique_50pct[[4]], H8_unique_50pct[[4]],
                              H9_unique_50pct[[4]], H10_unique_50pct[[4]], H11_unique_50pct[[4]], H12_unique_50pct[[4]],
                              H13_unique_50pct[[4]])
H_genes_gpt_50pct_R5 <- rbind(H1_unique_50pct[[5]], H2_unique_50pct[[5]], H3_unique_50pct[[5]], H4_unique_50pct[[5]],
                              H5_unique_50pct[[5]], H6_unique_50pct[[5]], H7_unique_50pct[[5]], H8_unique_50pct[[5]],
                              H9_unique_50pct[[5]], H10_unique_50pct[[5]], H11_unique_50pct[[5]], H12_unique_50pct[[5]],
                              H13_unique_50pct[[5]])


######## Each one of these individual control gene signatures will be used as input for AddmoduleScore and calculated similarly to the original gene signatures as shown in HallmarkScores.R script ####
