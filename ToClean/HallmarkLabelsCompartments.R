df_all$compartment <- cut(df_all$clusters, breaks = c(0,2.5,3.5,5), labels = c("Cancer", "Buffer", "TME"))


HallmarkActivity_labels <- df_all %>% select(sample, H1, H2, H3, H4, 
                  H5, H6,  H7, H8, H9, H10,  
                  H11, H12,  H13, compartment) %>%
  filter(compartment != "Buffer") %>% melt() %>%
  group_by(sample, compartment, variable) %>% 
  summarise(h_mean = mean(value), .groups = "drop") %>%
  group_by(sample, variable) %>% summarise(diff=h_mean[1]-h_mean[2]) %>%#diff cancer - TME 
  mutate(label=cut(diff, breaks = c(-Inf, -0.25, 0.25, Inf), labels = c("TME","Similar", "Cancer")))
  
HallmarkActivity_labels_TME <- filter(summs, variable %in% paste0("H", c(1,3,5,6,7,9,13)))
HallmarkActivity_labels_TME$label2 <- "Up" 
HallmarkActivity_labels_TME$label2[HallmarkActivity_labels_TME$label == "Cancer"] <- "Down"
HallmarkActivity_labels_TME$label2[HallmarkActivity_labels_TME$label == "Similar"] <- "Equal"
write.table(HallmarkActivity_labels_TME, "Desktop/HallmarkActivity_labels_TME.txt", sep = "\t")


HallmarkActivity_labels_Cancer <- filter(summs, variable %in% paste0("H", c(2,4,8,10,11,12)))
HallmarkActivity_labels_Cancer$label2 <- "Down" 
HallmarkActivity_labels_Cancer$label2[HallmarkActivity_labels_Cancer$label == "Cancer"] <- "Up"
HallmarkActivity_labels_Cancer$label2[HallmarkActivity_labels_Cancer$label == "Similar"] <- "Equal"
write.table(HallmarkActivity_labels_Cancer, "Desktop/HallmarkActivity_labels_Cancer.txt", sep = "\t")
