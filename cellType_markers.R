library(readr)
library(dplyr)
library(stringr)
#filter to get human specific genes
panglao <- read_tsv("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz")

## Immune system
panglao_immune <- panglao %>%  filter(str_detect(species,"Hs")
                                      & str_detect(`canonical marker`,"1")
                                      & str_detect(organ,"Immune system"))
panglao_immune <- panglao_immune %>%
  group_by(`cell type`) %>%
  summarise(geneset = list(`official gene symbol`))

panglao_immune <- setNames(panglao_immune$geneset, panglao_immune$`cell type`)

panglao_immune <- panglao_immune[sapply(panglao_immune, length) >= 10]


## Connective tissue
panglao_connective <- panglao %>%  filter(str_detect(species,"Hs")
                                          & str_detect(`canonical marker`,"1")
                                          & str_detect(organ,"Connective tissue"))
panglao_connective <- panglao_connective %>%
  group_by(`cell type`) %>%
  summarise(geneset = list(`official gene symbol`))

panglao_connective <- setNames(panglao_connective$geneset, panglao_connective$`cell type`)

panglao_connective <- panglao_connective[sapply(panglao_connective, length) >= 10]

## Vasculature tissue
panglao_Vasculature <- panglao %>%  filter(str_detect(species,"Hs")
                                           & str_detect(`canonical marker`,"1")
                                           & str_detect(organ,"Vasculature"))
panglao_Vasculature <- panglao_Vasculature %>%
  group_by(`cell type`) %>%
  summarise(geneset = list(`official gene symbol`))

panglao_Vasculature <- setNames(panglao_Vasculature$geneset, panglao_Vasculature$`cell type`)

panglao_Vasculature <- panglao_Vasculature[sapply(panglao_Vasculature, length) >= 10]

## Nervious system
panglao_nervous <- panglao %>%  filter(str_detect(species,"Hs")
                                           & str_detect(`canonical marker`,"1")
                                           & str_detect(organ,"Brain"))
panglao_nervous <- panglao_nervous %>%
  group_by(`cell type`) %>%
  summarise(geneset = list(`official gene symbol`))

panglao_nervous <- setNames(panglao_nervous$geneset, panglao_nervous$`cell type`)

panglao_nervous <- panglao_nervous[sapply(panglao_nervous, length) >= 10]

## Main TME
TME <- c(panglao_Vasculature, panglao_connective, panglao_immune, panglao_nervous)

saveRDS(TME, file = "~/rnaseq/expression/ST/Hallmarks/TME_cells.rds")
