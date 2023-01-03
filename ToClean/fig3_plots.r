
color_codes <-  list("#15CE59",  "#701717", "#CB3BBD", "#4969D4", "#E6880D", "#000000",  "#EE0F16", "#132892", "#8E909B", "#71189E","#05F3EB", "#890269", "#95641A")

# Rsquared increase vs Rsquared basal
lm_mat.basal.sub$diff <- lm_mat.1.sub$Rsquared - lm_mat.basal.sub$Rsquared
lm_mat.basal.sub$metric <-  abs(lm_mat.basal.sub$Rsquared*lm_mat.basal.sub$EstimateSlope)/lm_mat.basal.sub$EstimateSlope

ggplot(lm_mat.basal.sub, aes(x=diff, y=Rsquared)) +
  geom_point() + geom_text(data = lm_mat.basal.sub[lm_mat.basal.sub$diff > 0.05,], aes(label=meta), check_overlap = T) +
  theme_classic() + labs(x="Rsquared increase", y="Rsquared (basal)") + theme(axis.text = element_text(size = 10),
                                                                              axis.title = element_text(size=15))


ggplot(lm_mat.basal.sub[lm_mat.basal.sub$Ha==""], aes(x=Sample, y=Rsquared)) + geom_bar(stat="identity") + 
  theme_classic() + 
