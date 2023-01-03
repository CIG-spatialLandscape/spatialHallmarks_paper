


lm_mat <- data.frame(matrix(nrow = 0, ncol = 5))
for (sample in unique(mat$sample)) {
  print(sample)
  for (compartment in unique(mat$h_type)) {
    for (hallmark in unique(mat$hallmark)) {
      sub_mat <-  mat[mat$sample == sample & mat$hallmark == hallmark & mat$h_type == compartment, ]
      lm_out <- lm(h~d, sub_mat)
      lm_mat <- rbind(lm_mat, c(sample, hallmark, compartment, lm_out$coefficients[2],summary(lm_out)$adj.r.squared))
    }
  }
}


sub_lm <- lm_mat[lm_mat$Compartment=="TME" & lm_mat$Hallmark == "H3",]
colnames(sub_lm) <- c("Sample", "Hallmark", "Compartment", "Slope", "Rsquared")
sub_lm$Rsquared2 <-  as.numeric(sub_lm$Rsquared)*as.numeric(sub_lm$Slope)/abs(as.numeric(sub_lm$Slope))


tmp <- as.numeric(sub_lm$Rsquared2)
names(tmp) <- sub_lm$Sample
sub_lm$Sample <- factor(sub_lm$Sample, levels = names(sort(tmp)))
ggplot(sub_lm, aes(x=Sample, y=as.numeric(Rsquared2), fill=as.numeric(Rsquared2))) + geom_bar(stat="identity") + ylim(c(-0.3, 0.3))



############# Volcano

mat <- mat[mat$d <= 500, ]
lm_mat <- data.frame(matrix(nrow = 0, ncol = 6))
for (sample in unique(mat$sample)) {
  print(sample)
  for (compartment in unique(mat$h_type)) {
    for (hallmark in unique(mat$hallmark)) {
      sub_mat <-  mat[mat$sample == sample & mat$hallmark == hallmark & mat$h_type == compartment, ]
        if (nrow(sub_mat != 0)) {
        lm_out <- lm(h~d, sub_mat)
        lm_mat <- rbind(lm_mat, c(sample, hallmark, compartment, lm_out$coefficients[2],summary(lm_out)$adj.r.squared,summary(lm_out)$coefficients[2,4]))
        }
    }
  }
}

colnames(lm_mat) <- c("Sample", "Hallmark", "Compartment", "Slope", "Rsquared", "Pvalue")
lm_mat$Slope <- as.numeric(lm_mat$Slope)
lm_mat$Rsquared <- as.numeric(lm_mat$Rsquared)

ggplot(lm_mat[lm_mat$Hallmark=="H3",], aes(x=Slope, y=Rsquared, label=Sample)) + geom_point() + facet_wrap(~Compartment) + geom_text(aes(label=Sample))
ggplot(lm_mat, aes(x=Slope, y=Rsquared, col=Compartment)) + geom_point() + geom_text(aes(label=Sample))

#Parallel implementation
mat <- mat[mat$d <= 500, ]

library(parallel)
lm_mat <- data.frame(matrix(nrow = 0, ncol = 6))
for (sample in unique(mat$sample)) {
  print(sample)
  for (compartment in unique(mat$h_type)) {
    sub_mat <-  mat[mat$sample == sample & mat$h_type == compartment, ]
    cl <- makeCluster (13)  
    #copy the variables for each CPU
    clusterExport (cl, varlist = c("sub_mat", "sample", "compartment"))
    a <- parLapply(cl, paste0("H", 1:13), function(hallmark) {
      if (nrow(sub_mat != 0)) {
        lm_out <- lm(h~d, sub_mat[sub_mat$hallmark==hallmark,])
        return(c(sample, hallmark, compartment, lm_out$coefficients[2],summary(lm_out)$adj.r.squared,summary(lm_out)$coefficients[2,4]))
      }
    })
    stopCluster(cl)
    lm_mat <- rbind(lm_mat, do.call(rbind, a))
  }
}



files <- list.files("Desktop/IJC/datasets/IGTP/figuresPaper/distances/", full.names = T)
mat <- lapply(files, function(file){
  read.table(file, sep = "\t", header = T)
})
mat <- do.call(rbind, mat)

files <- sapply(list.files("Desktop/IJC/datasets/IGTP/figuresPaper/distances", full.names = F), function(x){
  strsplit(x, split = "_")[[1]][1]
})
library(parallel)


f <- function(sample) {
  print(sample)
  mat <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/distances/", sample, "_distances.txt"))
  mat <- mat[mat$d <= 500, ]
  lm_mat <- data.frame(matrix(nrow = 0, ncol = 6))
  for (compartment in unique(mat$h_type)) {
    for (hallmark in unique(mat$hallmark)) {
      sub_mat <-  mat[mat$sample == sample & mat$hallmark == hallmark & mat$h_type == compartment, ]
      if (nrow(sub_mat != 0)) {
        lm_out <- lm(h~d, sub_mat)
        lm_mat <- rbind(lm_mat, c(sample, hallmark, compartment, lm_out$coefficients[2],summary(lm_out)$adj.r.squared,summary(lm_out)$coefficients[2,4]))
      }
    }
  }
  return(lm_mat)
}
lm_list <- parLapply(cl, unique(files), f)

cl <- makeCluster (10) 
clusterExport (cl, varlist = c("files"))
#clusterExport (cl, varlist = c("mat"))
lm_list <- parLapply(cl, unique(files), function(sample) {
  message(sample)
  #print(sample)
  mat <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/distances/", sample, "_distances.txt"))
  mat <- mat[mat$d <= 500, ]
  lm_mat <- data.frame(matrix(nrow = 0, ncol = 6))
  for (compartment in unique(mat$h_type)) {
    for (hallmark in unique(mat$hallmark)) {
      sub_mat <-  mat[mat$sample == sample & mat$hallmark == hallmark & mat$h_type == compartment, ]
      if (nrow(sub_mat != 0)) {
        lm_out <- lm(h~d, sub_mat)
        lm_mat <- rbind(lm_mat, c(sample, hallmark, compartment, lm_out$coefficients[2],summary(lm_out)$adj.r.squared,summary(lm_out)$coefficients[2,4]))
      }
    }
  }
  return(lm_mat)
})

cl <- makeCluster (16)  
#copy the variables for each CPU
clusterExport (cl, varlist = c("counts", "l"))
#compute pairwise correlations between all spots
correlations <- parSapply(cl, 1:(l-1), function(i) {
  sapply((i+1):l, function(j) {
    cor(counts[,i], counts[,j])
  })
})
stopCluster(cl)

