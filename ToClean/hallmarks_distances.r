################################################################################

                        #Compute scores for each sample

################################################################################


start <- Sys.time() #compute execution time

#extract all the Seurat objects from the folder
files <- list.files(path="Desktop/enhanced", pattern="*.rds", full.names=TRUE, recursive=FALSE)

for (file in files) {
  print(file) #print in which step (file) we are
  STobject <- readRDS(file) #load RDS object file
  #extract the name of the sample, if the path is changed, the numbers have to be changed
  name <- strsplit(strsplit(file, "/")[[1]][3], "_")[[1]][1] 
  coord <- STobject@images[[1]]@coordinates[,3:4] #extract coordinates
  #compute distance matrix from coordinates (euclidean)
  dist_coord <- as.matrix(dist(coord, upper = T, diag = T)) 
  #set diagonal as 1 to avoid issues (dividing by 0)
  diag(dist_coord) <- 1 
  hallmarks <- STobject@meta.data[,paste0("H", 1:13)] #extract hallmarks scores
  #scale (and center) hallmarks data
  hallmarks <- scale(hallmarks) 
  #create empty matrix for the comparision scores
  scores <- matrix(0, nrow=13, ncol=13)
  rownames(scores) <- colnames(hallmarks) 
  colnames(scores) <- colnames(hallmarks)
  
  #iterate for each hallmark
  for (h1 in 1:ncol(hallmarks)) {
    #iterate for each hallmark to be compared (skip those already done)
    for (h2 in h1:ncol(hallmarks)) {
      
      x <- c()
      #set as reference each spot of the tissue 
      for (subspot in colnames(STobject)){
        #compute the sum of the values of dividing the hallmark 2 by the euclidean distance
        #and multiply the sum by the  score of hallmark 1 from the reference spot
        x <- c(x, hallmarks[subspot, h1]*mean(hallmarks[,h2]/dist_coord[rownames(hallmarks),subspot]))  
      }
      #store the mean of the reference spot values
      scores[h1, h2] <- mean(x)
    }
  }
  #save the matrix to a txt file
  write.table(scores, paste0("Desktop/matrix_1/", name, ".txt"), sep = "\t")
  #clean memory
  rm(STobject, coord, dist_coord, hallmarks)
  gc()
}
end <- Sys.time() #compute execution time
end-start #compute execution time


################################################################################

                    #Summarise and plot scores

################################################################################

#extract all matrix files from the path
files <- files <- list.files(path="Desktop/matrix_1", pattern="*.txt", full.names=TRUE, recursive=FALSE)
files <- files[-c(6, 7, 8, 10)]
#create a list of matrices
matrices <- list()
for (file in files) {
  #extract the name of the sample, if the path is changed, the numbers have to be changed
  name <- strsplit(strsplit(file, "/")[[1]][3], ".txt")[[1]][1]
  #load the matrix and add it to the list
  matrices[[name]] <- read.table(file, header = T, sep = "\t")
  #create symmetrix matrix (copy upper triangle to the lower traingle)
  matrices[[name]][lower.tri(matrices[[name]])] = t(matrices[[name]])[lower.tri(matrices[[name]])]
}

#scale each matrix separately
matrices <- lapply(matrices, scale)


#### Plot boxplots ####

#silly function
unify <- function(x) {
  x
}
#unify all matrices from the list
matrices_unified <- sapply(matrices, unify)

library(reshape2)
#create a long format table from the transpose matrix (comparisons as columns)
data <- melt(t(matrices_unified))

# add metadata to the table for plotting proposes 
h1 <- c()
h2 <- c()
for (h in colnames(matrices[[1]])) {
  h1 <- c(h1, rep(h, 8*13))
  h2 <- c(h2, rep(h, 8))
}
h2 <- rep(h2, 13)
data$h1 <- factor(h1,  paste0("H", 1:13))
data$h2 <- factor(h2,  paste0("H", 1:13))

#plot boxplot for each comparision
ggplot(data, aes(h2, value, fill=h2)) + geom_boxplot() + facet_grid(rows =  vars(h1)) +
  geom_hline(aes(yintercept=0), col="red", linetype = 'dotted') + theme_bw() + 
  theme(legend.position  = "none", axis.title = element_blank())



#### Plot heatmap ####
library(pheatmap)

#summarise by mean
tmp_all <- do.call(cbind, matrices) #combine all matrices from the list
dim(tmp_all) <- c(13,13,12) #set dimensions, last dimension correspond to number of samples used
mean_mat <- apply(tmp_all, c(1,2), mean) #summarise by mean

#set colnames and rownames to the hallmarks
colnames(mean_mat) <- paste0("H", 1:13)
rownames(mean_mat) <- paste0("H", 1:13)

pheatmap(mean_mat)


#summarise by median
tmp_all <- do.call(cbind, matrices) #combine all matrices from the list
dim(tmp_all) <- c(13,13,8) #set dimensions, last dimension correspond to number of samples used
median_mat <- apply(tmp_all, c(1,2), median) #summarise by meadian

#set colnames and rownames to the hallmarks
colnames(median_mat) <- hallmark_names
rownames(median_mat) <- hallmark_names

pheatmap(median_mat, cluster_rows = F, cluster_cols = T, clustering_distance_cols = "euclidean", angle_col = 315, border_color = NA)

#### Plot network ####
library(igraph)

#create a fully connected graph with 13 vertices (hallmarks)
fg <- make_full_graph(13)
fg <- set.vertex.attribute(fg, "name",  value=paste0("H", 1:13))

#extract weights of the edges from score matrix lower triangle
w <- c()
j <- 0
for (i in 1:12) {
  for (j in (i+1):13) {
    w <- c(w, median_mat[j, i]) #decide which method to summarize
  }
}

#add weights to the graph
E(fg)$weight <- w
#distribute vertices by the distance score (higher value, shorter)
coords <- layout_with_fr(fg, weights = w)

plot(fg, layout=coords, edge.width=w)







files <- files <- list.files(path="Desktop/matrix", pattern="*.txt", full.names=TRUE, recursive=FALSE)
files <- files[-c(6, 7, 8, 10)]
#create a list of matrices
matrices <- list()
for (file in files) {
  #extract the name of the sample, if the path is changed, the numbers have to be changed
  name <- strsplit(strsplit(file, "/")[[1]][3], ".txt")[[1]][1]
  #load the matrix and add it to the list
  matrices[[name]] <- read.table(file, header = T, sep = "\t")
  colnames(matrices[[name]]) <- hallmark_names
  rownames(matrices[[name]]) <- hallmark_names
  #create symmetrix matrix (copy upper triangle to the lower traingle)
  matrices[[name]][lower.tri(matrices[[name]])] = t(matrices[[name]])[lower.tri(matrices[[name]])]
  matrices[[name]] = matrices[[name]][c(2,4,8, 11, 12), -c(2,4,8,9,13,11, 12)]
  
}

#scale each matrix separately
matrices <- lapply(matrices, scale)







tmp_all <- do.call(cbind, matrices) #combine all matrices from the list
dim(tmp_all) <- c(5,6,8) #set dimensions, last dimension correspond to number of samples used
median_mat <- apply(tmp_all, c(1,2), median) #summarise by meadian

colnames(median_mat) <- hallmark_names[-c(2,4,8,9,13, 11, 12)]
rownames(median_mat) <- hallmark_names[c(2,4,8, 11, 12)]

pheatmap(median_mat, 
         cluster_rows = T, cluster_cols = T, 
         clustering_distance_cols = "euclidean",
         angle_col = 315, border_color = NA, fontsize = 20, color = viridis::inferno(100))


