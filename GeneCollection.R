##################################################
## Project: Cancer Hallmarks
  ## Script purpose: Obtain hallmark signatures using Pathway Commons
## Date: 22/12/2022
## Author: Sergi Cervilla & Mustafa Sibai
##################################################


library(tidyverse)
library(dplyr)
library(VennDiagram)

# New PathwayCommons database
######
#load pathwayCommons database
Paths <- read.delim("../PathwayCommons12.All.hgnc.txt")
# remove pairs of genes that does not belong to a pathway name
Paths <- Paths[-which(Paths$PATHWAY_NAMES == ""), ]
# split each entry of "PATHWAY_NAMES" column by ";" to obtain single pathway names
Paths <- separate_rows(Paths, 6, sep = ";")
# remove duplicate pairs and pathways
Paths <- Paths[duplicated(Paths) == F, ]
# transform into wide format by concatenating genes within a pathway
Paths_sum <- Paths %>%
  select(PARTICIPANT_A, PARTICIPANT_B, PATHWAY_NAMES) %>%
  group_by(PATHWAY_NAMES) %>%
  summarise(gene_list = unique(paste0(PARTICIPANT_A, ", ", PARTICIPANT_B, collapse = ", ")))

# remove CHEBI symbols (small molecules)
Paths_sum$gene_list <- sapply(Paths_sum$gene_list, function(x) paste(unique(unlist(str_split(x, ", "))[!grepl("^CHEBI:", unlist(str_split(x, ", ")))]), collapse = ", "))
# remove pathways that are empty
Paths_sum <- Paths_sum[-which(Paths_sum$gene_list == ""), ]

# count how many unique genes are found in each pathway
gene_lengths <- count.fields(textConnection(Paths_sum$gene_list), sep = ",")
Paths_sum$gene_lengths <- gene_lengths

# Filter by very low and very high number of genes per pathway
Paths_sum <- filter(Paths_sum, gene_lengths >= 4)
Paths_sum <- filter(Paths_sum, gene_lengths < 475)



## Serach for the keywords in each pathway to obtain hallmarks

# H1: Sustaining Proliferative Signaling
H1.kw <- c(
  "Cell Cycle", "Cell Division", "Growth Hormone", "EGF", "ERBB", "FGF", "PI3", "MAPK", "signaling by AKT1", "RAS", "MTOR", "MYC", "CDK2", "JAK/STAT", "STAT3",
  "signaling by FGFR1", "signaling by FGFR2", "signaling by FGFR3", "signaling by FGFR4"
)

# adding regular expression to the keywords
H1.kw <- paste0("\\b", H1.kw, "\\b")

# Collect the path containing the keywords
H1.pths <- Paths_sum %>%
  filter(grepl(paste(H1.kw, collapse = "|"), PATHWAY_NAMES, ignore.case = TRUE, ))

# concatenate all genes coming from H1 pathway
H1.genes <- data.frame(unlist(str_split(H1.pths$gene_list, ", ")))
H1.genes <- data.frame(H1.genes[duplicated(H1.genes) == FALSE, ])
colnames(H1.genes) <- "H1_genes"

H1.pths$hallmark <- "H1"
# check if BRAF is included in the gene collection
"BRAF" %in% H1.genes$H1_genes # TRUE

# H2: Evading growth supressors
H2.kw <- c("RB1", "P53", "APC", "SMAD4", "Cyclin", "Degradation of Cyclin", "Inactivation of Cyclin", "Checkpoint")

# adding regular expression to the keywords
H2.kw <- paste0("\\b", H2.kw, "\\b")

# Collect the path containing the keywords
H2.pths <- Paths_sum %>%
  filter(grepl(paste(H2.kw, collapse = "|"), PATHWAY_NAMES, ignore.case = TRUE))

# concatenate all genes coming from H2 pathway
H2.genes <- data.frame(unlist(str_split(H2.pths$gene_list, ", ")))
H2.genes <- data.frame(H2.genes[duplicated(H2.genes) == FALSE, ])
colnames(H2.genes) <- "H2_genes"

H2.pths$hallmark <- "H2"

# check if PTEN is included in the gene collection
"PTEN" %in% H2.genes$H2_genes # TRUE

# H3: Avoiding Immune Destruction
H3.kw <- c("immune", "IFN", "APC", "MHC", "PD-1 signaling")
# adding regular expression to the keywords
H3.kw <- paste0("\\b", H3.kw, "\\b")
# Collect the path containing the keywords
H3.pths <- Paths_sum %>%
  filter(grepl(paste(H3.kw, collapse = "|"), PATHWAY_NAMES, ignore.case = TRUE))
# concatenate all genes coming from H3 pathway
H3.genes <- data.frame(unlist(str_split(H3.pths$gene_list, ", ")))
H3.genes <- data.frame(H3.genes[duplicated(H3.genes) == FALSE, ])
colnames(H3.genes) <- "H3_genes"

# check if NLRC5 is included in the gene collection
"NLRC5" %in% H3.genes$H3_genes # FALSE


# Look for pathways that include NLRC5 in the gene list
H3.NLRC5 <- Paths_sum %>%
  filter(grepl(paste("\\<NLRC5,|\\<NLRC5[[:digit:]],|\\<NLRC5\\>|\\<NLRC5[[:digit:]]\\>", collapse = "|"), gene_list, ignore.case = TRUE))

H3.genes.NLRC5 <- data.frame(unlist(str_split(H3.NLRC5$gene_list, ", ")))
H3.genes.NLRC5 <- data.frame(H3.genes.NLRC5[duplicated(H3.genes.NLRC5) == FALSE, ])
colnames(H3.genes.NLRC5) <- "H3_genes_NLRC5"

# add NLRC5 to H3 gene list
colnames(H3.genes.NLRC5) <- colnames(H3.genes)
H3.genes <- rbind(H3.genes, H3.genes.NLRC5)

# combine paths
H3.pths <- rbind(H3.pths, H3.NLRC5)
H3.pths$hallmark <- "H3"



# H4: Enabling Replicative Immortality
H4.kw <- c("Telomere", "Recognition of DNA damage")
# adding regular expression to the keywords
H4.kw <- paste0("\\b", H4.kw, "\\b")

# Collect the paths containing the keywords
H4.pths <- Paths_sum %>%
  filter(grepl(paste(H4.kw, collapse = "|"), PATHWAY_NAMES, ignore.case = TRUE))
# concatenate all genes coming from H4 pathway
H4.genes <- data.frame(unlist(str_split(H4.pths$gene_list, ", ")))
H4.genes <- data.frame(H4.genes[duplicated(H4.genes) == FALSE, ])
colnames(H4.genes) <- "H4_genes"

H4.pths$hallmark <- "H4"

# H5: Tumor-promting inflammation
H5.kw <- c("inflammation", "SMAD", "chemokine", "cytokine", "interferon", "IFN")
# adding regular expression to the keywords
H5.kw <- paste0("\\b", H5.kw, "\\b")
# Collect the paths containing the keywords
H5.pths <- Paths_sum %>%
  filter(grepl(paste(H5.kw, collapse = "|"), PATHWAY_NAMES, ignore.case = TRUE))
# concatenate all genes coming from H5 pathway
H5.genes <- data.frame(unlist(str_split(H5.pths$gene_list, ", ")))
H5.genes <- data.frame(H5.genes[duplicated(H5.genes) == FALSE, ])
colnames(H5.genes) <- "H5_genes"

H5.pths$hallmark <- "H5"


# H6: Activating Invasion and Metastasis
# keywords with regular expression
H6.kw <- c("\\bEMT\\b", "E-Cadherin", "Metalloprot", "integrin", "\\bHedgehog signaling\\b")
# Collect the paths containing the keywords
H6.pths <- Paths_sum %>%
  filter(grepl(paste(H6.kw, collapse = "|"), PATHWAY_NAMES, ignore.case = TRUE))
# concatenate all genes coming from H6 pathway
H6.genes <- data.frame(unlist(str_split(H6.pths$gene_list, ", ")))
H6.genes <- data.frame(H6.genes[duplicated(H6.genes) == FALSE, ])
colnames(H6.genes) <- "H6_genes"

H6.pths$hallmark <- "H6"

# H7: Inducing Angiogenesis
# keywords with regular expression
H7.kw <- c("\\bWnt Signaling Pathway\\b", "\\bHypoxia\\b", "Angiogen", "\\bTie2\\b", "\\bPDGF\\b", "\\bNOTCH\\b", "\\bHIF\\b", "\\bVEGF\\b")
# Collect the paths containing the keywords
H7.pths <- Paths_sum %>%
  filter(grepl(paste(H7.kw, collapse = "|"), PATHWAY_NAMES, ignore.case = TRUE))
# concatenate all genes coming from H7 pathway
H7.genes <- data.frame(unlist(str_split(H7.pths$gene_list, ", ")))
H7.genes <- data.frame(H7.genes[duplicated(H7.genes) == FALSE, ])
colnames(H7.genes) <- "H7_genes"

H7.pths$hallmark <- "H7"

# H8: Genome Instability and mutation
H8.kw <- c("Dna repair", "DNA damage")
# adding regular expression to the keywords
H8.kw <- paste0("\\b", H8.kw, "\\b")
# Collect the paths containing the keywords
H8.pths <- Paths_sum %>%
  filter(grepl(paste(H8.kw, collapse = "|"), PATHWAY_NAMES, ignore.case = TRUE))
# concatenate all genes coming from H8 pathway
H8.genes <- data.frame(unlist(str_split(H8.pths$gene_list, ", ")))
H8.genes <- data.frame(H8.genes[duplicated(H8.genes) == FALSE, ])
colnames(H8.genes) <- "H8_genes"


"BRCA1" %in% H8.genes$H8_genes # TRUE

"BRCA2" %in% H8.genes$H8_genes # FALSE


# Look for pathways that include BRCA2 in the gene list
H8.BRCA <- Paths_sum %>%
  filter(grepl(paste("\\<BRCA2,|\\<BRCA2[[:digit:]],|\\<BRCA2\\>|\\<BRCA2[[:digit:]]\\>", collapse = "|"), gene_list, ignore.case = TRUE))

H8.genes.BRCA <- data.frame(unlist(str_split(H8.BRCA$gene_list, ", ")))
H8.genes.BRCA <- data.frame(H8.genes.BRCA[duplicated(H8.genes.BRCA) == FALSE, ])
colnames(H8.genes.BRCA) <- "H8_genes_BRCA"

table(H8.genes.BRCA$H8_genes_BRCA %in% H8.genes$H8_genes)

# add BRCA2 to H8 gene list
colnames(H8.genes.BRCA) <- colnames(H8.genes)
H8.genes <- rbind(H8.genes, H8.genes.BRCA)
H8.genes <- data.frame(H8.genes[duplicated(H8.genes) == FALSE, ])
colnames(H8.genes) <- "H8_genes"


# combine paths
H8.pths <- rbind(H8.pths, H8.BRCA)

H8.pths$hallmark <- "H8"


# H9: Resisting Cell Death
H9.kw <- c("Apoptosis", "Death", "BCL", "FAS", "TNF", "TRAIL", "CD95")
# adding regular expression to the keywords
H9.kw <- paste0("\\b", H9.kw, "\\b")
# Collect the paths containing the keywords
H9.pths <- Paths_sum %>%
  filter(grepl(paste(H9.kw, collapse = "|"), PATHWAY_NAMES, ignore.case = TRUE))
# concatenate all genes coming from H9 pathway
H9.genes <- data.frame(unlist(str_split(H9.pths$gene_list, ", ")))
H9.genes <- data.frame(H9.genes[duplicated(H9.genes) == FALSE, ])
colnames(H9.genes) <- "H9_genes"

H9.pths$hallmark <- "H9"

# H10: Deregulating cellular energetics
H10.kw <- c("HIF", "Glycolysis", "Insulin", "Warburg Effect", "AMPK", "Energy", "Mitochondrial Beta-Oxidation")
# adding regular expression to the keywords
H10.kw <- paste0("\\b", H10.kw, "\\b")
# Collect the paths containing the keywords
H10.pths <- Paths_sum %>%
  filter(grepl(paste(H10.kw, collapse = "|"), PATHWAY_NAMES, ignore.case = TRUE))
# concatenate all genes coming from H10 pathway
H10.genes <- data.frame(unlist(str_split(H10.pths$gene_list, ", ")))
H10.genes <- data.frame(H10.genes[duplicated(H10.genes) == FALSE, ])
colnames(H10.genes) <- "H10_genes"

H10.pths$hallmark <- "H10"

# H11: Senescent cells
H11.kw <- c("Senescence")
# Collect the paths containing the keywords
H11.pths <- Paths_sum %>%
  filter(grepl(paste(H11.kw, collapse = "|"), PATHWAY_NAMES, ignore.case = TRUE))
# Remove Senescent pathway that is involved in H4
H11.pths <- H11.pths %>%
  filter(!grepl(paste("Telomere", collapse = "|"), PATHWAY_NAMES, ignore.case = TRUE))
# concatenate all genes coming from H11 pathway
H11.genes <- data.frame(unlist(str_split(H11.pths$gene_list, ", ")))
H11.genes <- data.frame(H11.genes[duplicated(H11.genes) == FALSE, ])
colnames(H11.genes) <- "H11_genes"

H11.pths$hallmark <- "H11"

# H12: Nonmutational epigenetic reporgramming
# Regular expressions for the keywords
H12.kw <- c("\\bhistones?\\b", "\\bchromatin\\b")
# Collect the paths containing the keywords
H12.pths <- Paths_sum %>%
  filter(grepl(paste(H12.kw, collapse = "|"), PATHWAY_NAMES, ignore.case = TRUE))
# concatenate all genes coming from H12 pathway
H12.genes <- data.frame(unlist(str_split(H12.pths$gene_list, ", ")))
H12.genes <- data.frame(H12.genes[duplicated(H12.genes) == FALSE, ])
colnames(H12.genes) <- "H12_genes"

H12.pths$hallmark <- "H12"

# H13: Unlocking Phenotypic Plasticity
H13.kw <- c("stem cell", "pluripotent", "repress", "progenitor")
# adding regular expression to the keywords
H13.kw <- paste0("\\b", H13.kw, "\\b")
# Collect the paths containing the keywords
H13.pths <- Paths_sum %>%
  filter(grepl(paste(H13.kw, collapse = "|"), PATHWAY_NAMES, ignore.case = TRUE))
# concatenate all genes coming from H13 pathway
H13.genes <- data.frame(unlist(str_split(H13.pths$gene_list, ", ")))
H13.genes <- data.frame(H13.genes[duplicated(H13.genes) == FALSE, ])
colnames(H13.genes) <- "H13_genes"


H13.pths$hallmark <- "H13"


########## combine all hallmark genes ############

H1.genes$H <- "H1"
colnames(H1.genes)[1] <- "gene"
H2.genes$H <- "H2"
colnames(H2.genes)[1] <- "gene"
H3.genes$H <- "H3"
colnames(H3.genes)[1] <- "gene"
H4.genes$H <- "H4"
colnames(H4.genes)[1] <- "gene"
H5.genes$H <- "H5"
colnames(H5.genes)[1] <- "gene"
H6.genes$H <- "H6"
colnames(H6.genes)[1] <- "gene"
H7.genes$H <- "H7"
colnames(H7.genes)[1] <- "gene"
H8.genes$H <- "H8"
colnames(H8.genes)[1] <- "gene"
H9.genes$H <- "H9"
colnames(H9.genes)[1] <- "gene"
H10.genes$H <- "H10"
colnames(H10.genes)[1] <- "gene"
H11.genes$H <- "H11"
colnames(H11.genes)[1] <- "gene"
H12.genes$H <- "H12"
colnames(H12.genes)[1] <- "gene"
H13.genes$H <- "H13"
colnames(H13.genes)[1] <- "gene"

# final list of genes per hallmark
H.genes <- rbind(
  H1.genes, H2.genes, H3.genes, H4.genes, H5.genes, H6.genes, H7.genes, H8.genes,
  H9.genes, H10.genes, H11.genes, H12.genes, H13.genes
)

# number of unique genes
length(unique(H.genes$gene))
# number of genes per hallmark
table(H.genes$H)


# Concatenate pathways (final list of pathways per hallmark)
paths_final <- rbind(
  H1.pths,
  H2.pths,
  H3.pths,
  H4.pths,
  H5.pths,
  H6.pths,
  H7.pths,
  H8.pths,
  H9.pths,
  H10.pths,
  H11.pths,
  H12.pths,
  H13.pths
)
paths_final <- paths_final[4]

