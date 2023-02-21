library(STutility)
library(stringr)
library(dplyr)

path_objects <- "~/projects/Hallmarks/objects_mts/"

source("Desktop/IJC/git/spatialHallmarks_paper/utils/CoordinatesEnhanced.R") 
args = commandArgs(trailingOnly=TRUE)
sample <- args[1]
compartment <- args[2]

if (compartment == "Neoplastic") {clusters <- 1:2
} else if (compartment == "TME") {clusters <- 4:5
} else {clusters <- 1:5}



# original function: https://github.com/jbergenstrahle/STUtility/blob/master/R/Spatial_genes.R
#modification to apply it on enhance resolution
GetSpatNet_2 <- function (
  sample = NULL,
  spots = NULL,
  nNeighbours = NULL,
  maxdist = NULL,
  minK = 0
) {
  
  coordinates <- read.table(paste0(path_objects, sample, "_coords.txt"))
  coordinates$sample <- sample
  #compute real distances on ennhanced object
  coordinates <- subspot_coord(coordinates)
  
  
  # spatial information
  xys = setNames(coordinates[, c("realrow", "realcol", "sample")], c("x", "y", "sample"))
  xys <- xys[spots, ]
  
  
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
  hallmarks <- read.table(paste0("~/projects/Hallmarks/df_outs/", sample, ".txt"))
  data.use <- hallmarks %>% filter(clusters %in% cluster) %>% select(H1:H13) 
  data.use <- data.use[,features,drop=FALSE] %>% t()
  
  # Create a combined network for the samples
  CN <- do.call(rbind, GetSpatNet_2(sample = sample, spots = colnames(data.use), nNeighbours = nNeighbours, maxdist = maxdist))
  
  ###
  resCN <- as.matrix(data.frame(reshape2::dcast(CN, formula = from ~ to, value.var = "distance", fill = 0), row.names = 1))
  #resCN[resCN > 0] <- 1
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
  
  res <- data.frame(hallmark = rownames(data.use), cor = sp.cor, stringsAsFactors = F)
  res <- res[order(sp.cor, decreasing = T), ]
  rownames(res) <- res$hallmark
  return(res)
}

#compute Moran's for all hallmarks
df <- CorSpatialGenes_2(sample, cluster = clusters, features = paste0("H", 1:13), nNeighbours = 10)
write.table(df, paste0("~/projects/Hallmarks/SpAutocorrelation/Neoplastic/",sample, "_HallmarkMoran.txt"), sep = "\t")
