library(progeny)
library(tidyverse)

add_path_activities <- function(visium_slide, 
                                species = "human",
                                top = 1000, 
                                verbose = F,
                                assay = "SCT"){
  
  if(species == "mouse"){
    model <- progeny::getModel(organism = "Mouse", top = top)
    common_genes <- intersect(rownames(GetAssayData(visium_slide, assay = assay)), rownames(model))
    progeny_scores <- scale(t(progeny_scores))
    visium_slide[['progeny']] <- CreateAssayObject(counts = t(progeny_scores))
    
    #progeny_scores <- progeny::progeny(expr = GetAssayData(visium_slide, assay = assay),
    #                                  scale=TRUE, 
    #                                  organism="Mouse", 
    #                                  top=top, 
    #                                  perm=1, 
    #                                  verbose = verbose)
    
    #visium_slide[['progeny']] <- CreateAssayObject(counts = t(progeny_scores))
    
  }else if(species == "human"){
    
    model <- progeny::getModel(organism = "Human", top = top)
    common_genes <- intersect(rownames(GetAssayData(visium_slide, assay = assay)), rownames(model))
    progeny_scores <- t(model)[, common_genes] %*% GetAssayData(visium_slide, assay = assay)[common_genes, ]
    progeny_scores <- scale(t(progeny_scores))
    
    visium_slide[['progeny']] <- CreateAssayObject(counts = t(progeny_scores))
    
    #progeny_scores <- progeny::progeny(expr = GetAssayData(visium_slide, assay = assay),
    #                         scale=TRUE, 
    #                         organism="Human", 
    #                         top=top, 
    #                         perm=1,
    #                         verbose = verbose)
    
    #visium_slide[['progeny']] <- CreateAssayObject(data = t(progeny_scores))
    
  }
  
  return(visium_slide)
}

OC.st <- add_path_activities(OC.st)

pheatmap(mat = OC.st@assays$progeny@counts)
DoHeatmap(OC.st, features = rownames(OC.st))



progeny_scores_df <- 
  as.data.frame(t(GetAssayData(OC.st, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 

## We match Progeny scores with the cell clusters.
CellsClusters <- data.frame(Cell = names(Idents(OC.st)), 
                            CellType = as.character(Idents(OC.st)),
                            stringsAsFactors = FALSE)
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)


paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14, 
                        fontsize_row = 10, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "PROGENy (1000)", angle_col = 45,
                        treeheight_col = 20,  border_color = NA)



###3
gene_list <- read.table("Desktop/annotLookup.xls", header = T)
protein_genes <- gene_list$external_gene_name[gene_list$gene_biotype=="protein_coding"]

