#### Load data ####
files <- list.files("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/output_df", full.names = F)
files <- stringr::str_remove(files, pattern = ".txt")
df_all <- lapply(files, function(x) {
  df <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/output_df/", x, ".txt"), sep = "\t")
  hallmark <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/objects_mts/", x, "_hallmarks.txt"), sep = "\t")
  hallmark <- hallmark[-1,]
  df[,paste0("H", 1:13)] <- hallmark
  return(df)
})
df_all <- do.call(rbind, df_all)



df_all$spotid <- sapply(rownames(df_all), function(spot){
  x <- strsplit(spot,  ".", fixed = T)
  paste0(x[[1]][1], ".", substr(x[[1]][2],1, 1))
})


df_tme <- df_all[df_all$clusters %in% c(4,5),]
df_tme[, c("H2_cancer", "H4_cancer", "H8_cancer", "H10_cancer","H11_cancer", "H12_cancer")] <- NA
## Hallmark contribution 

for (x in unique(df_tme$sample)) {
  print(x)
  coord <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/objects_mts/", x, "_coords.txt"), sep = "\t")
  coord <- coord[df_all$spotid[df_all$sample==x],]
  #put real distances
  for (subspot in rownames(coord)) {
    realrow <-  coord[paste0(substr(subspot, 0, nchar(subspot)-2), ".5"), "row"]
    realrow <- (realrow+1)*100
    realcol <- trunc(coord[paste0(substr(subspot, 0, nchar(subspot)-2), ".3"), "col"])
    realcol <- (realcol+1)*100
    n <- substr(subspot, nchar(subspot), nchar(subspot))
    
    factor <- 55/3
    if (n == 1) {
      realrow <- realrow + factor
      realcol <- realcol + factor
    }
    else if (n == 2) {
      realrow <- realrow + factor
      realcol <- realcol - factor
    }
    else if (n == 3) {
      realrow <- realrow - factor
      realcol <- realcol + factor
    }
    else if (n == 4) {
      realrow <- realrow - factor
      realcol <- realcol - factor
    }
    else if (n == 5) {
      realrow <- realrow
      realcol <- realcol + 2 * factor
    }
    else if (n == 6) {
      realrow <- realrow
      realcol <- realcol - 2 * factor
    }
    coord[subspot, c("realrow", "realcol")] <- c(realrow, realcol)
    
    
  }
  
  tme_spots <- df_all$spotid[df_all$sample==x & df_all$clusters %in% c(4,5)]
  cancer_spots <- df_all$spotid[df_all$sample==x & df_all$clusters %in% c(1,2)]
  coord <- coord[c(tme_spots, cancer_spots), c("realrow", "realcol")]
  distances <- as.matrix(dist(coord))
  distances <- distances[tme_spots, cancer_spots]
  
  for (hallmark in c("H2", "H4", "H8", "H10","H11", "H12")) {
    print(hallmark)
    df_tme[df_tme$sample == x, paste0(hallmark, "_cancer")] <- sapply(tme_spots, function(spot){
      sum(1/distances[spot,]*df_all[df_all$sample==x & df_all$clusters %in% c(1,2), hallmark])
    })
  }
  rm(distances, coord)
  gc()
  
}

write.table(df_tme, "Desktop/df_tme58.txt", sep = "\t", quote = F)
gc()
.rs.restartR()

df_tme <- read.table("Desktop/df_tme.txt", sep = "\t")
df_tme2 <- read.table("Desktop/df_cancer.txt", sep = "\t")

df_tme2$spotid <- sapply(rownames(df_tme2), function(spot){
  x <- strsplit(spot,  ".", fixed = T)
  paste0(x[[1]][1], ".", substr(x[[1]][2],1, 1))
})

df_tme <- rbind(df_tme2, df_tme)

df_tme <- df_tme[df_tme$sample=="P7",]





print(hallmark)
hallmark <- "H2"

table(distances[tme_spots[127],] < 500)
neigh_spots <- cancer_spots[distances[tme_spots[127],] < 500]
length(1/distances[spot,]*df_all[df_all$sample==x & df_all$clusters %in% c(1,2), hallmark])

a <- sapply(tme_spots, function(spot){
  neigh_spots <- cancer_spots[distances[spot,] < 500]
  mean(1/distances[spot, neigh_spots]*df_all[df_all$sample==x & df_all$spotid %in% neigh_spots, hallmark])
})
