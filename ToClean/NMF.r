library(ica)
library(NNLM)
library(parallel)
FactorGeneLoadingPlot <- function (
  object,
  factor = 1,
  topn = 20,
  dark.theme = FALSE
) {
  
  # Check if NMF has been computed
  if (!"NMF" %in% names(object@reductions)) stop("NMF has not been computed ... \n", call. = FALSE)
  
  ftr <- paste0("factor_", factor)
  nmf <- object@reductions$NMF@feature.loadings[, ftr]
  gene <- names(nmf)
  df <- data.frame(gene, val = nmf, stringsAsFactors = F)
  df <- df[order(df$val, decreasing = T), ]
  df <- df[1:topn, ]
  df$gene <- factor(df$gene, levels = df$gene)
  p <- ggplot(df[1:topn, ], aes(reorder(gene, val), val)) +
    geom_bar(stat = "identity", fill = ifelse(dark.theme, "dark gray", "lightgray"), color = ifelse(dark.theme, "lightgray", "black"), width = 0.7) +
    coord_flip() +
    labs(x = "gene", y = "value")
  
  if (dark.theme) p <- p + DarkTheme() else p <- p + theme_minimal()
  
  return(p)
}


#' Run Non-negative Matrix Factorization
#'
#' Decompose an expression matrix A with non-negative elements into matrices WxH, also with
#' non-negative elements. W is the feature loading matrix (features x factors) and H is the
#' low dimensional embedding of the spots (factors x spots).
#'
#' @param object Seurat object
#' @param assay Assay Name of Assay NMF is being run on
#' @param features Features to compute the NMF for. Note that these features must be present in the
#' slot used to compute the NMF. By default, the `features` is set to `VariableFeatures(object)`
#' to include the most variable features selected in the normalization step.
#' @param nfactors Total Number of factors to compute and store (20 by default)
#' @param rescale Rescale data to make sure that values of the input matrix are non-n
#' @param reduction.name Dimensional reduction name, "NMF" by default
#' @param reduction.key Dimensional reduction key, specifies the prefix of the factor ids, e.g.
#' "factor_1", "factor_2", etc.
#' @param n.cores Number of threads to use in computation
#' @param order.by.spcor Order factors by spatial correlation
#'
#' @importFrom parallel detectCores
#'
#' @export
#'
RunNMF <- function (
  object,
  assay = NULL,
  slot = "scale.data",
  features = NULL,
  nfactors = 20,
  rescale = TRUE,
  reduction.name = "NMF",
  reduction.key = "factor_",
  n.cores = NULL,
  order.by.spcor = FALSE,
  sort.spcor.by.var = FALSE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  var.genes <- features %||% VariableFeatures(object)
  norm.counts <- GetAssayData(object, slot = slot, assay = assay)
  if (rescale) {
    norm.counts <- t(apply(norm.counts, 1, function(x) (x - min(x))/(max(x) - min(x))))
  }
  if (min(norm.counts) < 0) stop("Negative values are not allowed")
  nmf.results <- rnmf(A = norm.counts[var.genes, ], k = nfactors)
  #nmf.results$W <- swne::ProjectFeatures(norm.counts, nmf.results$H, n.cores = n.cores)
  feature.loadings <- nmf.results$W
  cell.embeddings <- t(nmf.results$H)
  
  # Set cores
  n.cores <- n.cores %||% {detectCores() - 1}
  
  # Order factors based on spatial correlation
  if (order.by.spcor) {
    CN <- do.call(rbind, GetSpatNet(object = object, nNeighbours = NULL, maxdist = NULL))
    resCN <- as.matrix(data.frame(reshape2::dcast(CN, formula = from ~ to, value.var = "distance", fill = 0), row.names = 1))
    resCN[resCN > 0] <- 1
    empty.CN <- matrix(0, nrow = nrow(cell.embeddings), ncol = nrow(cell.embeddings), dimnames = list(rownames(cell.embeddings), rownames(cell.embeddings)))
    colnames(resCN) <- gsub(pattern = "\\.", replacement = "-", x = colnames(resCN))
    colnames(resCN) <- gsub(pattern = "^X", replacement = "", x = colnames(resCN))
    empty.CN[rownames(resCN), colnames(resCN)] <- resCN
    listw <- mat2listw(empty.CN)
    fun <- function (x) lag.listw(listw, x, TRUE)
    
    # Calculate the lag matrix from the network
    tablag <- apply(cell.embeddings, 2, fun)
    
    # Split sp.cor by sample
    if (sort.spcor.by.var) {
      sp.cor.split <- do.call(rbind, lapply(unique(GetStaffli(object)@meta.data$sample), function(s) {
        tablag.split <- tablag[GetStaffli(object)@meta.data$sample == s, ]
        cell.embeddings.split <- cell.embeddings[GetStaffli(object)@meta.data$sample == s, ]
        unlist(lapply(1:ncol(cell.embeddings.split), function(i) {
          cor(tablag.split[, i], cell.embeddings.split[, i])
        }))
      }))
      order.vec <- order(apply(sp.cor.split, 2, var))
    } else {
      sp.cor <- unlist(lapply(1:ncol(cell.embeddings), function(i) {
        cor(cell.embeddings[, i], tablag[, i])
      }))
      order.vec <- order(sp.cor, decreasing = TRUE)
    }
    
    cell.embeddings <- cell.embeddings[, order.vec]
    colnames(cell.embeddings) <- paste0(reduction.key, 1:ncol(cell.embeddings))
  }
  
  rownames(x = feature.loadings) <- var.genes
  colnames(x = feature.loadings) <- paste0(reduction.key, 1:nfactors)
  rownames(x = cell.embeddings) <- colnames(x = object)
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
  reduction.data <- CreateDimReducObject (
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = assay,
    key = reduction.key
  )
  object[[reduction.name]] <- reduction.data
  return(object)
}



rnmf <- function (
  A,
  k,
  alpha = 0,
  init = "ica",
  n.cores = 1,
  loss = "mse",
  max.iter = 500,
  ica.fast = F
) {
  if (any(A < 0))
    stop("The input matrix contains negative elements !")
  if (k < 3)
    stop("k must be greater than or equal to 3 to create a viable SWNE plot")
  if (!init %in% c("ica", "nnsvd", "random")) {
    stop("Invalid initialization method")
  }
  A <- as.matrix(A)
  if (any(A < 0)) {
    stop("Input matrix has negative values")
  }
  if (init == "ica") {
    nmf.init <- ica_init(A, k, ica.fast = ica.fast)
  }
  else if (init == "nnsvd") {
    nmf.init <- nnsvd_init(A, k, LINPACK = T)
  }
  else {
    nmf.init <- NULL
  }
  if (is.null(nmf.init)) {
    nmf.res <- NNLM::nnmf(A, k = k, alpha = alpha, n.threads = n.cores,
                          loss = loss, max.iter = max.iter)
  }
  else {
    A.mean <- mean(A)
    zero.eps <- 1e-06
    nmf.init$W[nmf.init$W < zero.eps] <- 0
    nmf.init$H[nmf.init$H < zero.eps] <- 0
    zero.idx.w <- which(nmf.init$W == 0)
    zero.idx.h <- which(nmf.init$H == 0)
    nmf.init$W[zero.idx.w] <- runif(length(zero.idx.w), 0,
                                    A.mean/100)
    nmf.init$H[zero.idx.h] <- runif(length(zero.idx.h), 0,
                                    A.mean/100)
    nmf.res <- NNLM::nnmf(A, k = k, alpha = alpha, init = nmf.init,
                          n.threads = n.cores, loss = loss, max.iter = max.iter)
  }
  colnames(nmf.res$W) <- rownames(nmf.res$H) <- sapply(1:ncol(nmf.res$W),
                                                       function(i) paste("factor", i, sep = "_"))
  return(nmf.res)
}



#' Summarize features associated with cselected factors
#'
#' Extracts the top driving features per factor and returns
#'
#' @param object Seurat object
#' @param dims Factors to use
#' @param features.return Number of features to return per factor
#' @param features.use Select features (genes) to subset the data on
#'
#' @export
#'
SummarizeAssocFeatures <- function (
  object,
  dims = NULL,
  features.return = 10,
  features.use = NULL
) {
  
  if (!"NMF" %in% names(object@reductions)) stop(paste0("No factors available in Seurat object. Run RunNMF() first "), call. = FALSE)
  
  feature.factor.assoc <- object@reductions[["NMF"]]@feature.loadings
  if (!is.null(features.use)) {
    feature.factor.assoc <- feature.factor.assoc[features.use, ]
  }
  if (!is.null(dims)) {
    feature.factor.assoc <- feature.factor.assoc[, dims]
  }
  factor.features.df <- do.call("rbind", lapply(1:ncol(feature.factor.assoc),
                                                function(i) {
                                                  features.df <- data.frame(assoc_score = feature.factor.assoc[, i])
                                                  features.df$feature <- rownames(feature.factor.assoc)
                                                  features.df$factor <- colnames(feature.factor.assoc)[[i]]
                                                  features.df <- features.df[order(features.df$assoc_score, decreasing = T), ]
                                                  head(features.df, n = features.return)
                                                }))
  rownames(factor.features.df) <- NULL
  gene.loadings.selected <- feature.factor.assoc[unique(factor.features.df$feature), ]
  return(list(factor.features.df, gene.loadings.selected))
}


# Set a default value if an object is null
#
# @param lhs An object to set if it's null
# @param rhs The value to provide if x is null
#
# @return rhs if lhs is null, else lhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}
OC.st$
#' Initiate NMF using ICA
#'
#' @param A input matrix
#' @param k number of components to compute
#' @param ica.fast Should a fast implementation of ICA be used?
#'
#' @importFrom irlba irlba
#' @importFrom ica icafast

ica_init <- function (A, k, ica.fast = F)
{
  if (ica.fast) {
    pc.res.h <- irlba(t(A), nv = 50, maxit = 100,
                      center = rowMeans(A))
    ica.res.h <- icafast(pc.res.h$u, nc = k, maxit = 25,
                         tol = 1e-04)
    return(list(W = (A - Matrix::rowMeans(A)) %*% ica.res.h$S,
                H = t(ica.res.h$S)))
  }
  else {
    ica.res <- icafast(t(A), nc = k, maxit = 25, tol = 1e-04)
    return(list(W = ica.res$M, H = t(ica.res$S)))
  }
}

nnsvd_init <- function (A, k, LINPACK)
{
  size <- dim(A)
  m <- size[1]
  n <- size[2]
  W <- matrix(0, m, k)
  H <- matrix(0, k, n)
  s = svd(A, k, k, LINPACK = LINPACK)
  U <- s$u
  S <- s$d
  V <- s$v
  W[, 1] = sqrt(S[1]) * abs(U[, 1])
  H[1, ] = sqrt(S[1]) * abs(t(V[, 1]))
  for (i in seq(2, k)) {
    uu = U[, i]
    vv = V[, i]
    uup = .pos(uu)
    uun = .neg(uu)
    vvp = .pos(vv)
    vvn = .neg(vv)
    n_uup = .norm(uup)
    n_vvp = .norm(vvp)
    n_uun = .norm(uun)
    n_vvn = .norm(vvn)
    termp = n_uup %*% n_vvp
    termn = n_uun %*% n_vvn
    if (termp >= termn) {
      W[, i] = sqrt(S[i] * termp) * uup/n_uup
      H[i, ] = sqrt(S[i] * termp) * vvp/n_vvp
    }
    else {
      W[, i] = sqrt(S[i] * termn) * uun/n_uun
      H[i, ] = sqrt(S[i] * termn) * vvn/n_vvn
    }
  }
  return(list(W = W, H = H))
}


ST.DimPlot <- function (
  object,
  dims = c(1, 2),
  spots = NULL,
  indices = NULL,
  plot.type = "spots",
  blend = FALSE,
  min.cutoff = NA,
  max.cutoff = NA,
  pt.size = 1,
  pt.alpha = 1,
  pt.border = FALSE,
  reduction = NULL,
  palette = "MaYl",
  cols = NULL,
  dark.theme = FALSE,
  ncol = NULL,
  grid.ncol = NULL,
  center.zero = TRUE,
  channels.use = NULL,
  center.tissue = FALSE,
  verbose = FALSE,
  sb.size = 2.5,
  show.sb = TRUE,
  value.scale = c("samplewise", "all"),
  label.by = NULL,
  ...
) {
  
  # Check to see if Staffli object is present
  if (!"Staffli" %in% names(object@tools)) stop("Staffli object is missing from Seurat object. Cannot plot without coordinates", call. = FALSE)
  st.object <- object@tools$Staffli
  
  reduction <- reduction %||% {
    default.reductions <- c('umap', 'tsne', 'pca', 'ica')
    object.reductions <- FilterObjects(object = object, classes.keep = 'DimReduc')
    reduc.use <- min(which(x = default.reductions %in% object.reductions))
    default.reductions[reduc.use]
  }
  
  # prepare data
  signs <- sign(dims); dims <- abs(dims)
  spots <- spots %||% colnames(x = object)
  data <- Embeddings(object = object[[reduction]])[spots, dims, drop = FALSE]
  if (verbose) cat(paste0("Selected ", length(spots), " spots"))
  data <- as.data.frame(x = t(t(data)*signs))
  dims <- paste0(Key(object = object[[reduction]]), dims)
  
  # Select colorscale if palette is NULL
  palette.info <- palette.select(info = T)
  palette <- palette %||% {
    palette <- subset(palette.info, category == "div")$palette[1]
  }
  
  # Check that the number of dimensions are 2 or three if blending is active
  if (blend & !length(x = dims) %in% c(2, 3)) {
    stop(paste0("Blending dim plots only works with two or three dimensions. \n",
                "Number of dimensions provided: ", length(x = dims)), call. = F)
  }
  
  # Add group column to data
  data[,  "sample"] <- st.object[[spots, "sample", drop = TRUE]]
  
  # Obtain array coordinates
  image.type <- "empty"
  c(data, image.type) %<-% obtain.array.coords(st.object, data, image.type, spots)
  
  # Scale data values
  data <- feature.scaler(data, dims, min.cutoff, max.cutoff)
  
  # Subset by index
  if (!is.null(indices)) {
    if (!all(as.character(indices) %in% data[, "sample"])) stop(paste0("Index out of range. "), call. = FALSE)
    data <- data[data[, "sample"] %in% as.character(indices), ]
  } else {
    indices <- unique(data[,  "sample"]) %>% as.numeric()
  }
  
  # Fetch dims
  if (image.type != "empty") {
    dims.list <- lapply(st.object@dims, function(x) {x[2:3] %>% as.numeric()})
  } else {
    dims.list <- st.object@limits
  }
  
  # Subset dims by indices
  if (!is.null(indices)) dims.list <- dims.list[indices]
  
  # Prepare data for scalebar
  # --------------------------------------------------------------------
  pxum <- prep.sb(st.object, data, data.type, indices, FALSE, dims, dims.list, show.sb)
  # --------------------------------------------------------------------
  
  # Set feature scale limits
  value.scale <- match.arg(value.scale, c("samplewise", "all"))
  if (value.scale == "all") {
    limits <- c(min(data[, dims]), max(data[, dims]))
  } else if (value.scale == "samplewise") {
    limits <- NULL
  }
  
  # Add alternative label column
  if (!is.null(label.by)) {
    stopifnot(label.by %in% colnames(object@meta.data))
    stopifnot(class(object@meta.data[, label.by]) %in% c("character", "factor"))
    data <- setNames(cbind(data, object@meta.data[rownames(data), label.by]), nm = c(colnames(data), label.by))
  }
  
  # blend colors or plot each dimension separately
  if (blend) {
    colored.data <- apply(data[, 1:(ncol(data) - 3)], 2, scales::rescale)
    channels.use <- channels.use %||% c("red", "green", "blue")[1:ncol(colored.data)]
    
    if (verbose) cat(paste0("Blending colors for dimensions ",
                            paste0(ifelse(length(dims) == 2, paste0(dims[1], " and ", dims[2]), paste0(dims[1], dims[2], " and ", dims[2]))),
                            ": \n", paste(paste(dims, channels.use, sep = ":"), collapse = "\n")))
    
    spot.colors <- ColorBlender(colored.data, channels.use)
    data <- data[, (ncol(data) - 2):ncol(data)]
    p <- STPlot(data, data.type = "numeric", NULL, pt.size, pt.alpha, pt.border,
                palette, cols, ncol, spot.colors, center.zero, center.tissue,
                plot.title = paste(paste(dims, channels.use, sep = ":"), collapse = ", "),
                dims.list, split.labels = FALSE, pxum = pxum, sb.size = sb.size,
                dark.theme, NULL, label.by, ...)
    if (dark.theme) {
      p <- p + dark_theme()
    }
    return(p)
  } else {
    if (plot.type == "spots") {
      spot.colors <- NULL
      if (verbose) cat("Plotting dimensions:",
                       ifelse(length(dims) == 1, dims,  paste0(paste(dims[1:(length(dims) - 1)], collapse = ", "), " and ", dims[length(dims)])))
      
      plots <- lapply(X = dims, FUN = function(d) {
        ncol  <- ncol %||% 1
        plot <- STPlot(data, data.type = "numeric", d, pt.size, pt.alpha, pt.border,
                       palette, cols, ncol, spot.colors, center.zero, center.tissue,
                       d, dims.list, FALSE, pxum = pxum, sb.size = sb.size,
                       dark.theme, limits, label.by, ...)
        if (dark.theme) {
          plot <- plot + dark_theme()
        }
        return(plot)
      })
      
      # Draw plots
      grid.ncol <- grid.ncol %||% length(dims)
      p <- patchwork::wrap_plots(plots, ncol = grid.ncol)
      if (dark.theme) {
        p <- p & dark_theme()
      }
      return(p)
    } else if (plot.type == "smooth") {
      stop("Smooth options not yet implemented for dimred output ...")
      plots <- lapply(X = dims, FUN = function(d) {
        plot <- SmoothPlot(st.object, data, image.type, data.type = "numeric", d,
                           palette, cols, ncol, center.zero, dark.theme, highlight.edges, ...)
        return(plot)
      })
      
      # Draw plots
      grid.ncol <- grid.ncol %||% round(sqrt(length(x = plots)))
      grid.nrow <- ceiling(length(x = plots)/grid.ncol)
      
      stack <- c()
      for (i in 1:grid.nrow) {
        i <- i - 1
        stack <- c(stack, image_append(Reduce(c, plots[(i*grid.ncol + 1):(i*grid.ncol + grid.ncol)])))
      }
      
      final_img <- image_append(Reduce(c, stack), stack = T)
      print(final_img)
    }
  }
}

ica_init <- function(A, k, ica.fast = F) {
  if (ica.fast) {
    pc.res.h <- irlba::irlba(t(A), nv = 50, maxit = 500, center = rowMeans(A))
    ica.res.h <- ica::icafast(pc.res.h$u, nc = k, maxit = 50, tol = 1e-4)
    return(list(W = (A - Matrix::rowMeans(A)) %*% ica.res.h$S,
                H = t(ica.res.h$S)))
  } else {
    ica.res <- ica::icafast(t(A), nc = k, maxit = 50, tol = 1e-4)
    return(list(W = ica.res$M, H = t(ica.res$S)))
  }
}
