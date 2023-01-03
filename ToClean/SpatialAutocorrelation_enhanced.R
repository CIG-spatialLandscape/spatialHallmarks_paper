GetSpatNet_2 <- function (
  sample = NULL,
  spots = NULL,
  nNeighbours = NULL,
  maxdist = NULL,
  minK = 0
) {
  
  coordinates <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/objects_mts/",sample, "_coords.txt"))
  coordinates$sample <- sample
  coord <- coordinates
  ##### fix coordinates
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
    }  else if (n == 2) {
      realrow <- realrow + factor
      realcol <- realcol - factor
    }  else if (n == 3) {
      realrow <- realrow - factor
      realcol <- realcol + factor
    }  else if (n == 4) {
      realrow <- realrow - factor
      realcol <- realcol - factor
    }  else if (n == 5) {
      realrow <- realrow
      realcol <- realcol + 2 * factor
    }  else if (n == 6) {
      realrow <- realrow
      realcol <- realcol - 2 * factor
    }
    coord[subspot, c("realrow", "realcol")] <- c(realrow, realcol)
  }
  
 coordinates <- coord
  
  
  # spatial information
  xys = setNames(coordinates[, c("realrow", "realcol", "sample")], c("x", "y", "sample"))
  xys <- xys[spots, ]
  
  #xys <- xys[1:100,]
  
  # Split x, y, s
  xys.list <- split(xys, xys$sample)[unique(xys$sample)]
  
  # Obtain platforms
  platforms <- "Visium"
  
  # Compute network
  knn_spatial.norm.list <- lapply(seq_along(xys.list), function(i) {
    xys <- xys.list[[i]]
    
    # vector matching spot_ID and order
    spotnames <- rownames(xys)
    names(spotnames) <- c(1:nrow(xys)) %>% paste0()
    
    # Get spot distances
    sdist <- 100
    
    nNeighbours <- nNeighbours %||% ifelse(platforms[i] == "Visium", 6, 4)
    maxdist <- maxdist %||% ifelse(platforms[i] == "1k", 270*sdist, 150*sdist)
    
    if (!requireNamespace("dbscan")) install.packages("dbscan")
    knn_spatial <- dbscan::kNN(x = xys[, c("x", "y")] %>% as.matrix(), k = nNeighbours)
    knn_spatial.norm <- data.frame(from = rep(1:nrow(knn_spatial$id), nNeighbours),
                                   to = as.vector(knn_spatial$id),
                                   weight = 1/(1 + as.vector(knn_spatial$dist)),
                                   distance = as.vector(knn_spatial$dist))
    #nw_spatial.norm = igraph::graph_from_data_frame(knn_spatial.norm, directed = FALSE)
    #CN <- igraph::as_adjacency_matrix(nw_spatial.norm)
    
    # create network for coordinates
    spatnet <- knn_spatial.norm
    spatnet$from <- spotnames[spatnet$from]
    spatnet$to <- spotnames[spatnet$to]
    spatnet <- spatnet %>% group_by(from) %>% mutate(rnk = rank(distance)) %>% ungroup()
    spatnet =  subset(spatnet, distance <= maxdist | rnk <= minK)
    
    # Add coordinates
    spatnet <- cbind(spatnet, setNames(xys[spatnet$from, 1:2], paste0("start_", c("x", "y"))))
    spatnet <- cbind(spatnet, setNames(xys[spatnet$to, 1:2], paste0("end_", c("x", "y"))))
    
    return(spatnet)
  })
  
  return(knn_spatial.norm.list)
}


#' Find genes with high spatial autocorrelation
#'
#' This function can be used to find genes with spatial structure in ST datasets.
#' A more detailed decription of the algorithm is outlined in the Details section below.
#'
#' overview of method:
#' \itemize{
#'    \item{Build a connection network from the array x,y coordinates for each sample. For a 'Visium' array, this would typically be 6 neighbours
#'    because of the hexagonal structure of spots.}
#'    \item{Combine connection networks from multiple samples}
#'    \item{Compute the lag vector for each feature}
#'    \item{Compute the correlation between the lag vector and the original vector}
#' }
#' The connection network is build by defining edges between each spot and its `nNeighborurs` closest
#' neighbours that are within a maximum distance defined by `maxdist`. This is to make sure that spots
#' along the tissue edges or holes have the correct number of neighbours. A connection network is built for
#' each section separately but they are then combined into one large connection network so that the
#' autocorrelation can be computed for the whole dataset.
#'
#' Now that we have a neighbour group defined for
#' each spot, we can calculate the lag vector for each feature. The lag vector of a features is essentially
#' the summed expression of that feature in the neighbour groups, computed for all spots and can be thought
#' of as a "smoothing" estimate.
#'
#' If we consider a spot A and its neighbours nbA, a feature with high spatial
#' corelation should have similar expression levels in both groups. We can therefore compute the a
#' correlation score between the lag vector and the "normal" expression vector to get an estimate of
#' the spatial autocorrelation.
#'
#' @param object Seurat object
#' @param assay Name of assay the function is being run on
#' @param slot Slot to use as input [default: 'scale.data']
#' @param features Features to rank by spatial autocorrelation. If no features are provided, the
#' features are selected using the `VariableFeatures` function in Seurat, meaning that the top variable genes
#' will be used.
#' @param nNeighbours Number of neighbours to find for each spot, For Visium data, this parameter is set to
#' 6 because of the spots are arranged in a hexagonal pattern and should have maximum 6 neighbors.
#' @param maxdist Maximum allowed distance to define neighbouring spots [default: 1.5]. If not provided, a
#' maximum distance is automatically selected depending on the platform. For Visium data, this maximum distance
#' is set to 150 microns.
#'
#' @return data.frame with gene names and correlation scores
#'
#' @importFrom Matrix bdiag
#' @importFrom Seurat DefaultAssay GetAssayData VariableFeatures
#' @importFrom stats cor
#'
#' @export
#'
CorSpatialGenes_2 <- function (
  sample = NULL,
  cluster = NULL,
  slot = 'scale.data',
  features = NULL,
  nNeighbours = NULL,
  maxdist = NULL
) {
  
  if (!requireNamespace("spdep")) install.packages("spdep")
  
  
  
  # Obtain data
  hallmarks <- read.table(paste0("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/output_df/", sample, ".txt"))
  data.use <- hallmarks %>% filter(clusters %in% cluster) %>% select(H1:H13) %>% t()
  data.use <- data.use[features,]

  # Create a combined network for the samples
  CN <- do.call(rbind, GetSpatNet_2(sample = sample, spots = colnames(data.use), nNeighbours = nNeighbours, maxdist = maxdist))
  ###
  #data.use <- data.use[, rownames(coordinates[1:100,])]
  ###
  resCN <- as.matrix(data.frame(reshape2::dcast(CN, formula = from ~ to, value.var = "distance", fill = 0), row.names = 1))
  resCN[resCN > 0] <- 1
  empty.CN <- matrix(0, nrow = ncol(data.use), ncol = ncol(data.use), dimnames = list(colnames(data.use), colnames(data.use)))
  empty.CN[rownames(resCN), colnames(resCN)] <- resCN
  listw <- spdep::mat2listw(empty.CN)
  fun <- function (x) spdep::lag.listw(listw, x, TRUE)
  
  # Calculate the lag matrix from the network
  #tablag <- apply(t(data.use), 2, fun)
  tablag <- lapply(1:nrow(data.use), function(i) {
    fun(x = data.use[i, ])
  })
  tablag <- do.call(rbind, tablag)
  sp.cor <- unlist(lapply(1:nrow(data.use), function(i) {
    cor(data.use[i, ], tablag[i, ])
  }))
  
  res <- data.frame(gene = rownames(data.use), cor = sp.cor, stringsAsFactors = F)
  res <- res[order(sp.cor, decreasing = T), ]
  rownames(res) <- res$gene
  return(res)
}


df <- CorSpatialGenes_2("HCC2T", cluster = 4:5, features = paste0("H", c(1,3,5,6,7,9,13)))  %>% arrange(gene)
CorSpatialGenes_2("HCC2T", cluster = 1:2, features = paste0("H", c(2,4,8,10,11,12)))
rbind(data.frame(), df$gene)

tme <- data.frame()
samples <- list.files("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/output_df", full.names = F)
samples <- stringr::str_remove(samples, pattern = ".txt")
for (sample in samples) {
  print(sample)
  df <- CorSpatialGenes_2(sample, cluster = 4:5, features = paste0("H", c(1,3,5,6,7,9,13)), nNeighbours = 20, maxdist = 200)  %>% arrange(gene)
  tme <- rbind(tme, df$cor)
}
colnames(tme) <- df$gene
rownames(tme) <- samples
head(tme)


cancer <- data.frame()
samples <- list.files("Desktop/IJC/datasets/IGTP/figuresPaper/neighbours_experiment/output_df", full.names = F)
samples <- stringr::str_remove(samples, pattern = ".txt")
for (sample in samples) {
  print(sample)
  df <- CorSpatialGenes_2(sample, cluster = 1:2, features = paste0("H", c(2,4,8,10,11,12)), nNeighbours = 20, maxdist = 200)  %>% arrange(gene)
  cancer <- rbind(cancer, df$cor)
}
colnames(cancer) <- df$gene
rownames(cancer) <- samples
head(cancer)


tme.long <- tme %>% melt()
ggplot(tme.long, aes(x=variable, y=value, fill=variable)) + geom_boxplot()

cancer.long <- cancer %>% melt()
ggplot(cancer.long, aes(x=variable, y=value, fill=variable)) + geom_boxplot()

cancer.links <- read.table("Downloads/RF_direction_cancer.txt")
tme.links <- read.table("Downloads/RF_Imp_circos_direction.txt")

cancer_samples <- cancer.links  %>% filter(!is.na(link) & link != "")  %>% select(sample) %>% unique()
tme_samples <- tme.links  %>% filter(!is.na(link) & link != "")  %>% select(sample) %>% unique()

tme$link <- "None"
tme[tme_samples$sample, "link"] <- "LinkTME"
tme[cancer_samples$sample, "link"] <- "LinkCancer"
tme[intersect(cancer_samples$sample,tme_samples$sample), "link"] <- "<intersect"

tme.long <- tme %>% melt()
ggplot(tme.long, aes(x=link, y=value, fill=link)) + geom_boxplot() + facet_wrap(~variable)

cancer$link <- "None"
cancer[tme_samples$sample, "link"] <- "LinkTME"
cancer[cancer_samples$sample, "link"] <- "LinkCancer"
cancer[intersect(cancer_samples$sample,tme_samples$sample), "link"] <- "<intersect"

cancer.long <- cancer %>% melt()
ggplot(cancer.long, aes(x=link, y=value, fill=link)) + geom_boxplot() + facet_wrap(~variable)



CorSpatialGenes_2("HCC2T", cluster = 1:2, features = paste0("H", c(2,4,8,10,11,12)), nNeighbours = 10)
CorSpatialGenes_2("HCC2T", cluster = 1:2, features = paste0("H", c(2,4,8,10,11,12)), nNeighbours = 6)
CorSpatialGenes_2("HCC2T", cluster = 1:2, features = paste0("H", c(2,4,8,10,11,12)), nNeighbours = 16, maxdist = 200)
CorSpatialGenes_2("HCC2T", cluster = 1:2, features = paste0("H", c(2,4,8,10,11,12)), nNeighbours = 20, maxdist = 200)
CorSpatialGenes_2("HCC2T", cluster = 1:2, features = paste0("H", c(2,4,8,10,11,12)), nNeighbours = 25, maxdist = 200)

CorSpatialGenes_2("HCC2T", cluster = 1:2, features = paste0("H", c(2,4,8,10,11,12)), nNeighbours = 20, maxdist = 200)




df_tme <- read.table("Desktop/df_tme58.txt")

prox_tme <- df_tme %>% filter(sample == "HCC2T") %>% mutate(proximity_H4H1 = H1*H4_cancer, proximity_H4H6 = H6*H4_cancer) %>% 
  select(proximity_H4H6, proximity_H4H1) %>% melt()

ggplot(prox_tme, aes(x=variable, y=value)) + geom_boxplot()
plot(df_tme[df_tme$sample == "HCC2T", "H1"], df_tme[df_tme$sample == "HCC2T", "H4_cancer"])
plot(df_tme[df_tme$sample == "HCC2T", "H6"], df_tme[df_tme$sample == "HCC2T", "H4_cancer"])
cor(df_tme[df_tme$sample == "HCC2T", "H6"], df_tme[df_tme$sample == "HCC2T", "H4_cancer"])
cor(df_tme[df_tme$sample == "HCC2T", "H1"], df_tme[df_tme$sample == "HCC2T", "H4_cancer"])
