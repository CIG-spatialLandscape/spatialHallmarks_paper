RF_tme <- read.table("Downloads/RF_Imp_direction_TME.txt") %>% select(ID:predictor)
head(RF_tme)

RF_tme$links <- sapply(1:nrow(RF_tme), function(x) {
  ifelse(RF_tme[x, paste0(RF_tme$predictor[x], "_cor")] >= 0.6, "link", "None")
})
RF_tme$link_predictor <- "None"
RF_tme$link_predictor[RF_tme$links == "link"] <- paste0(RF_tme$predictor[RF_tme$links == "link"], "link")



RF_cancer <- read.table("Downloads/RF_Imp_direction_cancer.txt") %>% select(ID:predictor)
head(RF_cancer)

RF_cancer$links <- sapply(1:nrow(RF_cancer), function(x) {
  ifelse(RF_cancer[x, paste0(RF_cancer$predictor[x], "_cor")] >= 0.6, "link", "None")
})
RF_cancer$link_predictor <- "None"
RF_cancer$link_predictor[RF_cancer$links == "link"] <- paste0(RF_cancer$predictor[RF_cancer$links == "link"], "link")

# Check how many samples are for each response
RF_tme %>% filter(links == "link") %>% group_by(response) %>% summarise(n_samples = n())
RF_cancer %>% filter(links == "link") %>% group_by(response) %>% summarise(n_samples = n())

# Checl how many responses are for each sample
respxsample_tme <- RF_tme %>% filter(links == "link") %>% group_by(sample) %>% summarise(n_responses = n()) %>% arrange(-n_responses)
respxsample_cancer <- RF_cancer %>% filter(links == "link") %>% group_by(sample) %>% summarise(n_responses = n()) %>% arrange(-n_responses)

## Common most frequent samples
respxsample_tme.samples <- respxsample_tme %>% filter(n_responses >= quantile(n_responses, 0.75)) %>% select(sample)
respxsample_cancer.samples <- respxsample_cancer %>% filter(n_responses >= quantile(n_responses, 0.75)) %>% select(sample)

intersect(respxsample_tme.samples, respxsample_cancer.samples)

## Common least frequent samples
respxsample_tme.samples <- respxsample_tme %>% filter(n_responses <= quantile(n_responses, 0.25)) %>% select(sample)
respxsample_cancer.samples <- respxsample_cancer %>% filter(n_responses <= quantile(n_responses, 0.25)) %>% select(sample)

intersect(respxsample_tme.samples, respxsample_cancer.samples)

#
RF_all <- rbind(RF_tme %>% select(sample, type, response, predictor, links), RF_cancer %>% select(sample, type, response, predictor, links))
combinations <- combn(paste0("H", 1:13), 2)

combinations_count <- c()
for (i in 1:ncol(combinations)) {
  n <- 0
  for (sample in unique(RF_all$sample)) {
    predictors <- RF_all[RF_all$sample == sample & RF_all$links == "link", "predictor"]
    if (combinations[1,i] %in% predictors & combinations[2,i] %in% predictors) n <- n+1
  }
  combinations_count <- rbind(combinations_count, c(paste0(combinations[1,i], "_", combinations[2,i]), n))
}
combinations_count <- as.data.frame(combinations_count)
colnames(combinations_count) <- c("Hallmark pair", "n")
combinations_count$n <- as.numeric(combinations_count$n)

combinations_count$`Hallmark pair` <- factor(combinations_count$`Hallmark pair`, 
                                             levels = combinations_count %>% arrange(-n) %>% pull(`Hallmark pair`))

ggplot(combinations_count, aes(x=`Hallmark pair`, y=n)) + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 315, hjust = 0))

unique(iris$Species)
d <- iris[iris$Species != "Setosa",]
d$y <- "A"
d$y[d$Sepal.Width > mean(d$Sepal.Width)] <- "B"

createTable(compareGroups(y~Species, data=d, riskratio = T), show.ratio = T)



##
RF_all <- rbind(RF_tme %>% select(sample, type, response, predictor, links), RF_cancer %>% select(sample, type, response, predictor, links))
combinations <- combn(RF_tme %>% mutate(ResPre=paste0(response, "_", predictor)) %>% 
                        pull(ResPre) %>% unique(), 2)

combinations_count <- c()
for (i in 1:ncol(combinations)) {
  n <- 0
  for (sample in unique(RF_all$sample)) {
    predictors <- RF_all[RF_all$sample == sample & RF_all$links == "link", "predictor"]
    if (combinations[1,i] %in% predictors & combinations[2,i] %in% predictors) n <- n+1
  }
  combinations_count <- rbind(combinations_count, c(paste0(combinations[1,i], "_", combinations[2,i]), n))
}
combinations_count <- as.data.frame(combinations_count)
colnames(combinations_count) <- c("Hallmark pair", "n")
combinations_count$n <- as.numeric(combinations_count$n)

combinations_count$`Hallmark pair` <- factor(combinations_count$`Hallmark pair`, 
                                             levels = combinations_count %>% arrange(-n) %>% pull(`Hallmark pair`))

ggplot(combinations_count, aes(x=`Hallmark pair`, y=n)) + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 315, hjust = 0))

unique(iris$Species)
d <- iris[iris$Species != "Setosa",]
d$y <- "A"
d$y[d$Sepal.Width > mean(d$Sepal.Width)] <- "B"

createTable(compareGroups(y~Species, data=d, riskratio = T), show.ratio = T)


###

RF_tme %>% filter(links=="link") %>% mutate(ResPre=paste0(response, "_", predictor)) %>% group_by(type, ResPre) %>% count()
View(RF_tme %>% filter(links=="link") %>% mutate(ResPre=paste0(response, "_", predictor)) %>% group_by(ResPre) %>% count())

View(RF_tme %>% filter(links=="link") %>% mutate(ResPre=paste0(response, "_", predictor)) %>% group_by(ResPre) %>% summarise(n=n(), samples=paste0(sample, collapse = "; ")))


View(RF_cancer %>% filter(links=="link") %>% mutate(ResPre=paste0(response, "_", predictor)) %>% group_by(type, ResPre) %>% count())
View(RF_cancer %>% filter(links=="link") %>% mutate(ResPre=paste0(response, "_", predictor)) %>% group_by(ResPre) %>% summarise(n=n(), samples=paste0(sample, collapse = "; ")))

RF_cancer %>% filter(links=="link") %>% mutate(ResPre=paste0(response, "_", predictor)) %>% group_by(ResPre) %>% summarise(n=n(), samples=paste0(sample, collapse = "; ")) %>% slice_max(n, n=1) %>% pull(samples) %>% str_split(pattern = ";", simplify = T) %>% as.character()

intersect(RF_cancer %>% filter(links=="link") %>% mutate(ResPre=paste0(response, "_", predictor)) %>% group_by(ResPre) %>% summarise(n=n(), samples=paste0(sample, collapse = "; ")) %>% slice_max(n, n=1) %>% pull(samples) %>% str_split(pattern = ";", simplify = T) %>% as.character(),
          RF_tme %>% filter(links=="link") %>% mutate(ResPre=paste0(response, "_", predictor)) %>% group_by(ResPre) %>% summarise(n=n(), samples=paste0(sample, collapse = "; ")) %>% slice_max(n, n=1) %>% pull(samples) %>% str_split(pattern = ";", simplify = T) %>% as.character()
)



RF_all_ResPre_pos <- rbind(RF_tme %>% filter(links=="link") %>% mutate(ResPre=paste0(response, "_", predictor)) %>% select(sample, type, predictor, response, ResPre) %>% mutate(experiment="TME"),
                       RF_cancer %>% filter(links=="link") %>% mutate(ResPre=paste0(predictor, "_", response)) %>% select(sample, type, predictor, response, ResPre) %>% mutate(experiment="Cancer"))
View(RF_all_ResPre %>%  group_by(ResPre) %>% summarise(n=n(), TME=sum(experiment=="TME"), Cancer=sum(experiment=="Cancer")))

df <- RF_all_ResPre_pos %>%  group_by(ResPre) %>% summarise(n=n(), TME=sum(experiment=="TME"), Cancer=sum(experiment=="Cancer")) %>% select(-n) %>% melt()     
df$ResPre <- factor(df$ResPre, 
                    levels = RF_all_ResPre_pos %>%  group_by(ResPre) %>% summarise(n=n(), TME=sum(experiment=="TME"), Cancer=sum(experiment=="Cancer")) %>% arrange(-n) %>% pull(ResPre))

ggplot(df, aes(x=ResPre, y=value, fill=variable)) + geom_bar(stat = "identity") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 315, hjust = 0))



######### Negative

RF_tme_neg <- read.table("Downloads/RF_Imp_direction_TME.txt") %>% select(ID:predictor)
head(RF_tme_neg)

RF_tme_neg$links <- sapply(1:nrow(RF_tme_neg), function(x) {
  ifelse(RF_tme_neg[x, paste0(RF_tme_neg$predictor[x], "_cor")] <= -0.6, "link", "None")
})
RF_tme_neg$link_predictor <- "None"
RF_tme_neg$link_predictor[RF_tme_neg$links == "link"] <- paste0(RF_tme_neg$predictor[RF_tme_neg$links == "link"], "link")



RF_cancer_neg <- read.table("Downloads/RF_Imp_direction_cancer.txt") %>% select(ID:predictor)
head(RF_cancer)

RF_cancer_neg$links <- sapply(1:nrow(RF_cancer_neg), function(x) {
  ifelse(RF_cancer_neg[x, paste0(RF_cancer_neg$predictor[x], "_cor")] <= -0.6, "link", "None")
})
RF_cancer_neg$link_predictor <- "None"
RF_cancer_neg$link_predictor[RF_cancer_neg$links == "link"] <- paste0(RF_cancer_neg$predictor[RF_cancer_neg$links == "link"], "link")


RF_all_ResPre_neg <- rbind(RF_tme_neg %>% filter(links=="link") %>% mutate(ResPre=paste0(response, "_", predictor)) %>% select(sample, type, predictor, response, ResPre) %>% mutate(experiment="TME"),
                       RF_cancer_neg %>% filter(links=="link") %>% mutate(ResPre=paste0(predictor, "_", response)) %>% select(sample, type, predictor, response, ResPre) %>% mutate(experiment="Cancer"))

df <- RF_all_ResPre %>%  group_by(ResPre) %>% summarise(n=n(), TME=sum(experiment=="TME"), Cancer=sum(experiment=="Cancer")) %>% select(-n) %>% melt()     
df$ResPre <- factor(df$ResPre, 
                    levels = RF_all_ResPre %>%  group_by(ResPre) %>% summarise(n=n(), TME=sum(experiment=="TME"), Cancer=sum(experiment=="Cancer")) %>% arrange(-n) %>% pull(ResPre))

ggplot(df, aes(x=ResPre, y=value, fill=variable)) + geom_bar(stat = "identity") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 315, hjust = 0))




#####


pos_samples <- RF_all_ResPre_pos %>%  group_by(ResPre) %>% summarise(n=n(), TME=sum(experiment=="TME"), Cancer=sum(experiment=="Cancer"), samples = paste0(sample, collapse = ";")) %>% 
  arrange(-n) %>% slice_head(prop = 0.25) %>% pull(samples) %>% str_split(pattern = ";") %>% do.call(what = c)
neg_samples <- RF_all_ResPre_neg %>%  group_by(ResPre) %>% summarise(n=n(), TME=sum(experiment=="TME"), Cancer=sum(experiment=="Cancer"), samples = paste0(sample, collapse = ";")) %>% 
  arrange(-n) %>% slice_head(prop = 0.25) %>% pull(samples) %>% str_split(pattern = ";") %>% do.call(what = c)
intersect(pos_samples, neg_samples) #23



pos_samples_tme <- RF_tme %>% filter(links=="link") %>% mutate(ResPre=paste0(response, "_", predictor)) %>% select(sample, type, predictor, response, ResPre) %>%  group_by(ResPre) %>% summarise(n=n(), samples = paste0(sample, collapse = ";")) %>% 
  arrange(-n) %>% slice_head(prop = 0.25) %>% pull(samples) %>% str_split(pattern = ";") %>% do.call(what = c)
length(unique(pos_samples_tme))
pos_samples_cancer <- RF_cancer %>% filter(links=="link") %>% mutate(ResPre=paste0(response, "_", predictor)) %>% select(sample, type, predictor, response, ResPre) %>%  group_by(ResPre) %>% summarise(n=n(), samples = paste0(sample, collapse = ";")) %>% 
  arrange(-n) %>% slice_head(prop = 0.25) %>% pull(samples) %>% str_split(pattern = ";") %>% do.call(what = c)
length(unique(pos_samples_cancer))
intersect(pos_samples_cancer, pos_samples_tme) 
intersect(intersect(pos_samples, neg_samples), pos_samples_tme) #16
intersect(intersect(pos_samples, neg_samples), pos_samples_cancer) #9



neg_samples_tme <- RF_tme_neg %>% filter(links=="link") %>% mutate(ResPre=paste0(response, "_", predictor)) %>% select(sample, type, predictor, response, ResPre) %>%  group_by(ResPre) %>% summarise(n=n(), samples = paste0(sample, collapse = ";")) %>% 
  arrange(-n) %>% slice_head(prop = 0.25) %>% pull(samples) %>% str_split(pattern = ";") %>% do.call(what = c)
length(unique(neg_samples_tme))
neg_samples_cancer <- RF_cancer_neg %>% filter(links=="link") %>% mutate(ResPre=paste0(response, "_", predictor)) %>% select(sample, type, predictor, response, ResPre) %>%  group_by(ResPre) %>% summarise(n=n(), samples = paste0(sample, collapse = ";")) %>% 
  arrange(-n) %>% slice_head(prop = 0.25) %>% pull(samples) %>% str_split(pattern = ";") %>% do.call(what = c)
length(unique(neg_samples_cancer))
intersect(neg_samples_cancer, neg_samples_tme)
intersect(intersect(pos_samples, neg_samples), neg_samples_tme) #3
intersect(intersect(pos_samples, neg_samples), neg_samples_cancer) #11

##### Non-linear

RF_tme_nonlinear <- read.table("Downloads/RF_Imp_direction_TME.txt") %>% select(ID:predictor)
head(RF_tme_nonlinear)

RF_tme_nonlinear$links <- sapply(1:nrow(RF_tme_nonlinear), function(x) {
  ifelse(RF_tme_nonlinear[x, paste0(RF_tme_nonlinear$predictor[x], "_cor")] > -0.6 | RF_tme_nonlinear[x, paste0(RF_tme_nonlinear$predictor[x], "_cor")] < 0.6  , "link", "None")
})
RF_tme_nonlinear$link_predictor <- "None"
RF_tme_nonlinear$link_predictor[RF_tme_nonlinear$links == "link"] <- paste0(RF_tme_nonlinear$predictor[RF_tme_nonlinear$links == "link"], "link")



RF_cancer_nonlinear <- read.table("Downloads/RF_Imp_direction_cancer.txt") %>% select(ID:predictor)


RF_cancer_nonlinear$links <- sapply(1:nrow(RF_cancer_nonlinear), function(x) {
  ifelse(RF_cancer_nonlinear[x, paste0(RF_cancer_nonlinear$predictor[x], "_cor")] > -0.6 | RF_cancer_nonlinear[x, paste0(RF_cancer_nonlinear$predictor[x], "_cor")] < 0.6, "link", "None")
})
RF_cancer_nonlinear$link_predictor <- "None"
RF_cancer_nonlinear$link_predictor[RF_cancer_nonlinear$links == "link"] <- paste0(RF_cancer_nonlinear$predictor[RF_cancer_nonlinear$links == "link"], "link")


RF_all_ResPre_nonlinear <- rbind(RF_tme_nonlinear %>% filter(links=="link") %>% mutate(ResPre=paste0(response, "_", predictor)) %>% select(sample, type, predictor, response, ResPre) %>% mutate(experiment="TME"),
                           RF_cancer_nonlinear %>% filter(links=="link") %>% mutate(ResPre=paste0(predictor, "_", response)) %>% select(sample, type, predictor, response, ResPre) %>% mutate(experiment="Cancer"))

df <- RF_all_ResPre_nonlinear %>%  group_by(ResPre) %>% summarise(n=n(), TME=sum(experiment=="TME"), Cancer=sum(experiment=="Cancer")) %>% select(-n) %>% melt()     
df$ResPre <- factor(df$ResPre, 
                    levels = RF_all_ResPre_nonlinear %>%  group_by(ResPre) %>% summarise(n=n(), TME=sum(experiment=="TME"), Cancer=sum(experiment=="Cancer")) %>% arrange(-n) %>% pull(ResPre))

ggplot(df, aes(x=ResPre, y=value, fill=variable)) + geom_bar(stat = "identity") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 315, hjust = 0))




####

pos_samples <- RF_all_ResPre_pos %>%  group_by(ResPre) %>% summarise(n=n(), TME=sum(experiment=="TME"), Cancer=sum(experiment=="Cancer"), samples = paste0(sample, collapse = ";")) %>% 
  arrange(-n) %>% slice_head(prop = 0.25) %>% pull(samples) %>% str_split(pattern = ";") %>% do.call(what = c)
neg_samples <- RF_all_ResPre_neg %>%  group_by(ResPre) %>% summarise(n=n(), TME=sum(experiment=="TME"), Cancer=sum(experiment=="Cancer"), samples = paste0(sample, collapse = ";")) %>% 
  arrange(-n) %>% slice_head(prop = 0.25) %>% pull(samples) %>% str_split(pattern = ";") %>% do.call(what = c)
nonlinear_samples <- RF_all_ResPre_nonlinear %>%  group_by(ResPre) %>% summarise(n=n(), TME=sum(experiment=="TME"), Cancer=sum(experiment=="Cancer"), samples = paste0(sample, collapse = ";")) %>% 
  arrange(-n) %>% slice_head(prop = 0.25) %>% pull(samples) %>% str_split(pattern = ";") %>% do.call(what = c)

intersect(intersect(pos_samples, neg_samples), nonlinear_samples)

library(VennDiagram)
dev.off()
draw.triple.venn(area1 = length(unique(pos_samples)), area2 = length(unique(neg_samples)), 
                 area3 = length(unique(nonlinear_samples)), n12 = length(intersect(pos_samples, neg_samples)), 
                 n23 = length(intersect(nonlinear_samples, neg_samples)), n13 = length(intersect(pos_samples, nonlinear_samples)), 
                 n123 = length(intersect(intersect(pos_samples, neg_samples), nonlinear_samples)), category = c("Positive", "Negative", "Nonlinear"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))


num_samples_type <- RF_tme %>% select(sample, type) %>% unique() %>% group_by(type) %>% summarise(n=n())

inter <- intersect(intersect(pos_samples, neg_samples), nonlinear_samples)
inter_type <- sapply(inter, function(x){return(annotation_tissue[[x]])},USE.NAMES = F)
tbl1 <- round(table(inter_type) / filter(num_samples_type, type %in% sort(unique(inter_type))) %>% pull(n) * 100, 2)

inter_neg_non <- unique(neg_samples)[!unique(neg_samples) %in% inter]
inter_type_neg_non <- sapply(inter_neg_non, function(x){return(annotation_tissue[[x]])},USE.NAMES = F)
tbl2 <- round(table(inter_type_neg_non) / filter(num_samples_type, type %in% sort(unique(inter_type_neg_non))) %>% pull(n) * 100, 2)

inter_pos_non <- unique(pos_samples)[!unique(pos_samples) %in% inter]
inter_type_pos_non <- sapply(inter_pos_non, function(x){return(annotation_tissue[[x]])},USE.NAMES = F)
tbl3 <- round(table(inter_type_pos_non) / filter(num_samples_type, type %in% sort(unique(inter_type_pos_non))) %>% pull(n) * 100, 2)

inter_non <- unique(nonlinear_samples)[!unique(nonlinear_samples) %in% c(inter, inter_pos_non, inter_neg_non)]
inter_type_non <- sapply(inter_non, function(x){return(annotation_tissue[[x]])},USE.NAMES = F)
tbl4 <- round(table(inter_type_non) / filter(num_samples_type, type %in% sort(unique(inter_type_non))) %>% pull(n) * 100, 2)


col_tissue_map <- setNames(Seurat::DiscretePalette(15, palette = "polychrome") [c(1:3, 5, 6:11)], sort(unique(RF_tme$type)))
ggplot(as.data.frame(tbl1), aes(x=inter_type, y=Freq, fill=inter_type)) + geom_bar(stat = "identity") + theme_classic() +
  scale_fill_manual(values = col_tissue_map) + labs(x="") + theme(axis.text.x = element_text(angle = 315, hjust = 0),
                                                                  legend.position = "none") + ylim(c(0,100))
ggsave("Desktop/tbl1.png", bg = "white")
ggplot(as.data.frame(tbl2), aes(x=inter_type_neg_non, y=Freq, fill=inter_type_neg_non)) + geom_bar(stat = "identity") + theme_classic() +
  scale_fill_manual(values = col_tissue_map) + labs(x="") + theme(axis.text.x = element_text(angle = 315, hjust = 0),
                                                                  legend.position = "none") + ylim(c(0,100))
ggsave("Desktop/tbl2.png", bg = "white")
ggplot(as.data.frame(tbl3), aes(x=inter_type_pos_non, y=Freq, fill=inter_type_pos_non)) + geom_bar(stat = "identity") + theme_classic() +
  scale_fill_manual(values = col_tissue_map) + labs(x="") + theme(axis.text.x = element_text(angle = 315, hjust = 0),
                                                                  legend.position = "none") + ylim(c(0,100))
ggsave("Desktop/tbl3.png", bg = "white")
ggplot(as.data.frame(tbl4), aes(x=inter_type_non, y=Freq, fill=inter_type_non)) + geom_bar(stat = "identity") + theme_classic() +
  scale_fill_manual(values = col_tissue_map) + labs(x="") + theme(axis.text.x = element_text(angle = 315, hjust = 0),
                                                                  legend.position = "none") + ylim(c(0,100))
ggsave("Desktop/tbl4.png", bg = "white")
