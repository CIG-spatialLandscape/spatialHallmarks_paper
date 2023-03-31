##################################################
## Project: Cancer Hallmarks
## Script purpose: Modifications of BayesSpace and Seurat plotting to combine and plot high resolution sub-spot plots
## Author: Sergi Cervilla* & Mustafa Sibai*
##################################################

library(assertthat)
clusterPlot <- function(sce, label="spatial.cluster",
                        palette=NULL, color=NULL,
                        platform=NULL, is.enhanced=NULL,
                        ...) {
  
  if (is.null(platform))
    platform <- .bsData(sce, "platform", "Visium")
  if (is.null(is.enhanced))
    is.enhanced <- .bsData(sce, "is.enhanced", FALSE)
  
  vertices <- .make_vertices(sce, label, platform, is.enhanced)
  
  ## No borders around subspots by default
  if (is.null(color)) {
    color <- ifelse(is.enhanced, NA, "#d8dcd6")
  }
  
  splot <- ggplot(data=vertices, 
                  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~factor(fill))) +
    geom_polygon(color=color, ...) +
    labs(fill="Cluster") +
    coord_equal() +
    theme_void()
  
  if (!is.null(palette))
    splot <- splot + scale_fill_manual(values=palette)
  
  splot
}

#' Plot spatial gene expression.
#' 
#' @param sce SingleCellExperiment. If \code{feature} is specified and is a 
#'   string, it must exist as a row in the specified assay of \code{sce}.
#' @param feature Feature vector used to color each spot. May be the name of a
#'   gene/row in an assay of \code{sce}, or a vector of continuous values.
#' @param assay.type String indicating which assay in \code{sce} the expression
#'   vector should be taken from.
#' @param low,mid,high Optional hex codes for low, mid, and high values of the
#'   color gradient used for continuous spot values.
#' @param diverging If true, use a diverging color gradient in
#'   \code{featurePlot()} (e.g. when plotting a fold change) instead of a
#'   sequential gradient (e.g. when plotting expression).
#' @inheritParams spatialPlot
#' 
#' @return Returns a ggplot object.
#' 
#' @examples
#' sce <- exampleSCE()
#' featurePlot(sce, "gene_2")
#' 
#' @family spatial plotting functions
#'
#' @importFrom ggplot2 ggplot aes_ geom_polygon scale_fill_gradient scale_fill_gradient2 coord_equal labs theme_void
#' @importFrom scales muted
#' @importFrom assertthat assert_that
#' @export
featurePlot <- function(sce, feature,
                        assay.type="logcounts", 
                        diverging=FALSE,
                        low=NULL, high=NULL, mid=NULL,
                        color=NULL,
                        platform=NULL, is.enhanced=NULL,
                        ...) {
  
  if (is.null(platform))
    platform <- .bsData(sce, "platform", "Visium")
  if (is.null(is.enhanced))
    is.enhanced <- .bsData(sce, "is.enhanced", FALSE)
  
  ## extract expression from logcounts if a gene name is passed.
  ## otherwise, assume a vector of counts was passed and let
  ## .make_vertices helpers check validity
  if (is.character(feature)) {
    assert_that(feature %in% rownames(sce),
                msg=sprintf("Feature %s not in SCE.", feature))
    fill <- assay(sce, assay.type)[feature, ]
    fill.name <- feature
  } else {
    fill <- feature
    
    ## this could be an argument, but it's easily overwritten with labs()
    ## and we should encourage composing ggplot functions instead
    fill.name <- "Expression"
  }
  
  vertices <- .make_vertices(sce, fill, platform, is.enhanced)
  
  ## No borders around subspots by default
  if (is.null(color)) {
    color <- ifelse(is.enhanced, NA, "#d8dcd6")
  }
  
  splot <- ggplot(data=vertices, 
                  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +
    geom_polygon(color=color, ...) +
    labs(fill=fill.name) +
    coord_equal() +
    theme_void()
  
  if (diverging) {
    low <- ifelse(is.null(low), "#F0F0F0", low)
    high <- ifelse(is.null(high), muted("red"), high)
    splot <- splot + scale_fill_gradient(low=low, high=high)
  } else {
    low <- ifelse(is.null(low), muted("blue"), low)
    mid <- ifelse(is.null(mid), "#F0F0F0", mid)
    high <- ifelse(is.null(high), muted("red"), high)
    splot <- splot + scale_fill_gradient2(low=low, mid=mid, high=high)
  }
  
  splot
}



#' Make vertices outlining spots/subspots for geom_polygon()
#' 
#' @param sce SingleCellExperiment with row/col in colData
#' @param fill Name of a column in \code{colData(sce)} or a vector of values to
#'   use as fill for each spot
#' @param platform "Visium" or "ST", used to determine spot layout
#' @param is.enhanced If true, \code{sce} contains enhanced subspot data instead
#'   of spot-level expression. Used to determine spot layout.
#'   
#' @return Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
#'   vertices outlining the spot's border
#' 
#' @keywords internal
.make_vertices <- function(sce, fill, platform, is.enhanced) {
  cdata <- data.frame(colData(sce))
  
  if (platform == "Visium") {
    if (is.enhanced) {
      vertices <- .make_triangle_subspots(cdata, fill)
    } else {
      vertices <- .make_hex_spots(cdata, fill)
    }
  } else if (platform == "ST") {
    if (is.enhanced) {
      vertices <- .make_square_spots(cdata, fill, scale.factor=(1/3))
    } else {
      vertices <- .make_square_spots(cdata, fill)
    }
  } else {
    stop("Unsupported platform: \"", platform, "\". Cannot create spot layout.")
  }
  
  vertices
}

#' Helper to extract x, y, fill ID from colData
#' 
#' @return Dataframe of (x.pos, y.pos, fill) for each spot
#' 
#' @keywords internal
#' @importFrom assertthat assert_that
.select_spot_positions <- function(cdata, x="col", y="row", fill="spatial.cluster") {
  ## Provide either a column name or vector of labels/values
  assert_that(is.vector(fill) || is.character(fill) || is.factor(fill))
  
  ## I think this is the best way to check if something is a string
  if (is.character(fill) && length(fill) == 1) {
    spot_positions <- cdata[, c(x, y, fill)]
    colnames(spot_positions) <- c("x.pos", "y.pos", "fill")    
  } else if (is.vector(fill) || is.factor(fill)) {
    assert_that(nrow(cdata) == length(fill))
    spot_positions <- cdata[, c(x, y)]
    colnames(spot_positions) <- c("x.pos", "y.pos")    
    spot_positions$fill <- fill
  }
  spot_positions$spot <- rownames(spot_positions)
  
  spot_positions
}

#' Compute vertex coordinates for each spot in frame of plot
#'
#' @param spot_positions Center for hex, top left for square
#' @param vertex_offsets Data frame of (x, y) offsets wrt spot position for each
#'   vertex of spot
#' 
#' @return Cartesian product of positions and offsets, with coordinates
#'   computed as (pos + offset)
#'
#' @keywords internal
.make_spot_vertices <- function(spot_positions, vertex_offsets) {
  spot_vertices <- merge(spot_positions, vertex_offsets)
  spot_vertices$x.vertex <- spot_vertices$x.pos + spot_vertices$x.offset
  spot_vertices$y.vertex <- spot_vertices$y.pos + spot_vertices$y.offset
  
  as.data.frame(spot_vertices)
}

#' Make vertices for each hex spot
#' 
#' @return Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
#'   vertices outlining the spot's border
#' 
#' @keywords internal
.make_hex_spots <- function(cdata, fill) {
  ## R = circumradius, distance from center to vertex
  ## r = inradius, distance from center to edge midpoint
  r <- 1/2
  R <- (2 / sqrt(3)) * r
  
  spot_positions <- .select_spot_positions(cdata, fill=fill)
  spot_positions <- .adjust_hex_centers(spot_positions)
  
  ## vertices of each hex (with respect to center coordinates)
  ## start at top center, loop clockwise
  vertex_offsets <- data.frame(x.offset=c(0, r, r, 0, -r, -r),
                               y.offset=c(-R, -R/2, R/2, R, R/2, -R/2))
  
  spot_vertices <- .make_spot_vertices(spot_positions, vertex_offsets)
  
  ## Flip to match image orientation
  spot_vertices$y.vertex <- -spot_vertices$y.vertex
  
  spot_vertices
}

#' Adjust hex spot positions so hexagons are adjacent to each other in plot
#'
#' Spots are regular hexagons with one unit of horizontal distance
#' between centers
#' 
#' @return Shifted spot centers
#' 
#' @keywords internal
.adjust_hex_centers <- function(spot_positions) {
  ## R = circumradius, distance from center to vertex
  ## r = inradius, distance from center to edge midpoint
  r <- 1/2
  R <- (2 / sqrt(3)) * r
  
  ## Start at (1-indexed origin)
  spot_positions$x.pos <- spot_positions$x.pos - min(spot_positions$x.pos) + 1
  spot_positions$y.pos <- spot_positions$y.pos - min(spot_positions$y.pos) + 1
  
  ## Shift centers up so rows are adjacent
  spot_positions$y.pos <- spot_positions$y.pos * R * (3/2)
  
  ## Spot columns are offset by row
  ## (i.e. odd rows have odd numbered columns, even rows have even)
  ## Shift centers to the left so columns are adjacent (but hexes stay offset)
  spot_positions$x.pos <- (spot_positions$x.pos + 1) / 2
  
  spot_positions
}

#' Make vertices for each square spot
#'
#' Squares are simple, just mae a unit square at each array coordinate
#' 
#' @return Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
#'   vertices outlining the spot's border
#' 
#' @keywords internal
.make_square_spots <- function(cdata, fill="spatial.cluster", scale.factor=1) {
  spot_positions <- .select_spot_positions(cdata, fill=fill)
  
  vertex_offsets <- data.frame(x.offset=c(0, 1, 1, 0),
                               y.offset=c(0, 0, 1, 1))
  vertex_offsets <- vertex_offsets * scale.factor
  
  .make_spot_vertices(spot_positions, vertex_offsets)
}

#' Helper to pull out subspot position columns
#' Probably redundant with select_spot_positions above, but we need subspot.idx
#' 
#' @return Dataframe of (x.pos, y.pos, fill) for each spot
#' 
#' @keywords internal
.select_subspot_positions <- function(cdata, x="spot.col", y="spot.row", fill="spatial.cluster") {
  ## Provide either a column name or vector of labels/values
  assert_that(is.vector(fill) || is.character(fill) || is.factor(fill))
  
  if (is.character(fill) && length(fill) == 1) {
    spot_positions <- cdata[, c(x, y, "subspot.idx", fill)]
    colnames(spot_positions) <- c("x.pos", "y.pos", "subspot.idx", "fill")
  } else if (is.vector(fill) || is.factor(fill)) {
    assert_that(nrow(cdata) == length(fill))
    spot_positions <- cdata[, c(x, y, "subspot.idx")]
    colnames(spot_positions) <- c("x.pos", "y.pos", "subspot.idx")    
    spot_positions$fill <- fill
  }
  
  spot_positions$spot <- rownames(spot_positions)
  
  spot_positions
}

#' Make vertices for each triangle subspot of a hex
#' 
#' @return Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
#'   vertices outlining the spot's border
#'
#' @keywords internal
.make_triangle_subspots <- function(cdata, fill="spatial.cluster") {
  spot_positions <- .select_subspot_positions(cdata, x="spot.col", y="spot.row", fill=fill)
  spot_positions <- .adjust_hex_centers(spot_positions)
  
  ## R = circumradius, distance from center to vertex
  ## r = inradius, distance from center to edge midpoint
  r <- 1/2
  R <- (2 / sqrt(3)) * r
  
  ## Make lists of triangle vertices (with respect to hex center)
  ## subspot.idx is same ordering as `shift` in spatialEnhance
  ## that is, beginning in top right and proceeding clockwise, (1, 5, 3, 4, 6, 2)
  ## NOTE: however, we need to reflect over x-axis to match global negation of y-coordinate
  vertex_offsets <- do.call(rbind, list(
    data.frame(x.offset=c(0, 0, r), y.offset=c(0, -R, -R/2), subspot.idx=3),
    data.frame(x.offset=c(0, r, r), y.offset=c(0, -R/2, R/2), subspot.idx=5),
    data.frame(x.offset=c(0, r, 0), y.offset=c(0, R/2, R), subspot.idx=1),
    data.frame(x.offset=c(0, 0, -r), y.offset=c(0, R, R/2), subspot.idx=2),
    data.frame(x.offset=c(0, -r, -r), y.offset=c(0, R/2, -R/2), subspot.idx=6),
    data.frame(x.offset=c(0, -r, 0), y.offset=c(0, -R/2, -R), subspot.idx=4)
  ))
  
  ## note that instead of cartesian product, `merge()` does an outer join
  ## on subspot.idx here
  spot_vertices <- .make_spot_vertices(spot_positions, vertex_offsets)
  spot_vertices$y.vertex <- -spot_vertices$y.vertex
  
  spot_vertices
}

library(scales)
library(grid)
library(patchwork)
library(RColorBrewer)
library(Seurat)
library(ggplot2)
GeomSpatial <- ggproto(
  "GeomSpatial",
  Geom,
  required_aes = c("x", "y"),
  extra_params = c("na.rm", "image", "image.alpha", "crop"),
  default_aes = aes(
    shape = 21,
    colour = NA,
    point.size.factor = 1.0,
    fill = NA,
    alpha = NA,
    stroke = 0
  ),
  setup_data = function(self, data, params) {
    data <- ggproto_parent(Geom, self)$setup_data(data, params)
    # We need to flip the image as the Y coordinates are reversed
    data$y = max(data$y) - data$y + min(data$y)
    data
  },
  draw_key = draw_key_point,
  draw_panel = function(data, panel_scales, coord, image, image.alpha, crop) {
    # This should be in native units, where
    # Locations and sizes are relative to the x- and yscales for the current viewport.
    if (!crop) {
      y.transform <- c(0, nrow(x = image)) - panel_scales$y.range
      data$y <- data$y + sum(y.transform)
      panel_scales$x$continuous_range <- c(0, ncol(x = image))
      panel_scales$y$continuous_range <- c(0, nrow(x = image))
      panel_scales$y.range <- c(0, nrow(x = image))
      panel_scales$x.range <- c(0, ncol(x = image))
    }
    z <- coord$transform(
      data.frame(x = c(0, ncol(x = image)), y = c(0, nrow(x = image))),
      panel_scales
    )
    # Flip Y axis for image
    z$y <- -rev(z$y) + 1
    wdth <- z$x[2] - z$x[1]
    hgth <- z$y[2] - z$y[1]
    vp <- viewport(
      x = unit(x = z$x[1], units = "npc"),
      y = unit(x = z$y[1], units = "npc"),
      width = unit(x = wdth, units = "npc"),
      height = unit(x = hgth, units = "npc"),
      just = c("left", "bottom")
    )
    img.grob <- GetImage(object = image)
    
    img <- editGrob(grob = img.grob, vp = vp)
    # spot.size <- slot(object = image, name = "spot.radius")
    spot.size <- Radius(object = image)
    coords <- coord$transform(data, panel_scales)
    pts <- pointsGrob(
      x = coords$x,
      y = coords$y,
      pch = data$shape,
      size = unit(spot.size, "npc") * data$point.size.factor,
      gp = gpar(
        col = alpha(colour = coords$colour, alpha = coords$alpha),
        fill = alpha(colour = coords$fill, alpha = coords$alpha),
        lwd = coords$stroke)
    )
    vp <- viewport()
    gt <- gTree(vp = vp)
    if (image.alpha > 0) {
      if (image.alpha != 1) {
        img$raster = as.raster(
          x = matrix(
            data = alpha(colour = img$raster, alpha = image.alpha),
            nrow = nrow(x = img$raster),
            ncol = ncol(x = img$raster),
            byrow = TRUE)
        )
      }
      gt <- addGrob(gTree = gt, child = img)
    }
    gt <- addGrob(gTree = gt, child = pts)
    # Replacement for ggname
    gt$name <- grobName(grob = gt, prefix = 'geom_spatial')
    return(gt)
    # ggplot2:::ggname("geom_spatial", gt)
  }
)


geom_spatial <-  function(
  mapping = NULL,
  data = NULL,
  image = image,
  image.alpha = image.alpha,
  crop = crop,
  stat = "identity",
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  ...
) {
  layer(
    geom = GeomSpatial,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, image = image, image.alpha = image.alpha, crop = crop, ...)
  )
}


SingleSpatialPlot <- function (data, image, cols = NULL, image.alpha = 1, pt.alpha = NULL, 
                               crop = TRUE, pt.size.factor = NULL, stroke = 0.25, col.by = NULL, 
                               alpha.by = NULL, cells.highlight = NULL, cols.highlight = c("#DE2D26", 
                                                                                           "grey50"), geom = c("spatial", "interactive", "poly"), 
                               na.value = "grey50") 
{
  geom <- match.arg(arg = geom)
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find '", col.by, "' in data, not coloring", 
            call. = FALSE, immediate. = TRUE)
    col.by <- NULL
  }
  col.by <- col.by %iff% paste0("`", col.by, "`")
  alpha.by <- alpha.by %iff% paste0("`", alpha.by, "`")
  if (!is.null(x = cells.highlight)) {
    highlight.info <- SetHighlight(cells.highlight = cells.highlight, 
                                   cells.all = rownames(x = data), sizes.highlight = pt.size.factor, 
                                   cols.highlight = cols.highlight[1], col.base = cols.highlight[2])
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- "highlight"
    levels(x = data$ident) <- c(order, setdiff(x = levels(x = data$ident), 
                                               y = order))
    data <- data[order(data$ident), ]
  }
  plot <- ggplot(data = data, aes_string(x = colnames(x = data)[2], 
                                         y = colnames(x = data)[1], alpha = alpha.by))
  plot <- switch(EXPR = geom, spatial = {
    if (is.null(x = pt.alpha)) {
      plot <- plot + geom_spatial(point.size.factor = pt.size.factor, 
                                  data = data, image = image, image.alpha = image.alpha, 
                                  crop = crop, stroke = stroke, )
    } else {
      plot <- plot + geom_spatial(point.size.factor = pt.size.factor, 
                                  data = data, image = image, image.alpha = image.alpha, 
                                  crop = crop, stroke = stroke, alpha = pt.alpha)
    }
    plot + coord_fixed() + theme(aspect.ratio = 1)
  }, interactive = {
    plot + geom_spatial_interactive(data = tibble(grob = list(GetImage(object = image, 
                                                                       mode = "grob"))), mapping = aes_string(grob = "grob"), 
                                    x = 0.5, y = 0.5) + geom_point(mapping = aes_string(color = col.by)) + 
      xlim(0, ncol(x = image)) + ylim(nrow(x = image), 
                                      0) + coord_cartesian(expand = FALSE)
  }, poly = {
    data$cell <- rownames(x = data)
    data[, c("x", "y")] <- NULL
    data <- merge(x = data, y = GetTissueCoordinates(object = image, 
                                                     qhulls = TRUE), by = "cell")
    plot + geom_polygon(data = data, mapping = aes_string( 
      group = "cell")) + coord_fixed() + theme_cowplot()
  }, stop("Unknown geom, choose from 'spatial' or 'interactive'", 
          call. = FALSE))
  if (!is.null(x = cells.highlight)) {
    plot <- plot + scale_fill_manual(values = cols.highlight)
  }
  if (!is.null(x = cols) && is.null(x = cells.highlight)) {
    if (length(x = cols) == 1 && (is.numeric(x = cols) || 
                                  cols %in% rownames(x = brewer.pal.info))) {
      scale <- scale_fill_brewer(palette = cols, na.value = na.value)
    }
    else if (length(x = cols) == 1 && (cols %in% c("alphabet", 
                                                   "alphabet2", "glasbey", "polychrome", "stepped"))) {
      colors <- DiscretePalette(length(unique(data[[col.by]])), 
                                palette = cols)
      scale <- scale_fill_manual(values = colors, na.value = na.value)
    }
    else {
      scale <- scale_fill_manual(values = cols, na.value = na.value)
    }
    #plot <- plot + scale
  }
  plot <- plot + NoAxes() + theme(panel.background = element_blank())
  return(plot)
}

SpatialPlot <- function (object, group.by = NULL, features = NULL, images = NULL, 
                         cols = NULL, image.alpha = 1, crop = TRUE, slot = "data", 
                         min.cutoff = NA, max.cutoff = NA, cells.highlight = NULL, 
                         cols.highlight = c("#DE2D26", "grey50"), facet.highlight = FALSE, 
                         label = FALSE, label.size = 5, label.color = "white", label.box = TRUE, 
                         repel = FALSE, ncol = NULL, combine = TRUE, pt.size.factor = 1.6, 
                         alpha = c(1, 1), stroke = 0.25, interactive = FALSE, do.identify = FALSE, 
                         identify.ident = NULL, do.hover = FALSE, information = NULL) 
{
  if (isTRUE(x = do.hover) || isTRUE(x = do.identify)) {
    warning("'do.hover' and 'do.identify' are deprecated as we are removing plotly-based interactive graphics, use 'interactive' instead for Shiny-based interactivity", 
            call. = FALSE, immediate. = TRUE)
    interactive <- TRUE
  }
  if (!is.null(x = group.by) & !is.null(x = features)) {
    stop("Please specific either group.by or features, not both.")
  }
  images <- images %||% Images(object = object, assay = DefaultAssay(object = object))
  if (length(x = images) == 0) {
    images <- Images(object = object)
  }
  if (length(x = images) < 1) {
    stop("Could not find any spatial image information")
  }
  if (is.null(x = features)) {
    if (interactive) {
      return(ISpatialDimPlot(object = object, image = images[1], 
                             group.by = group.by, alpha = alpha))
    }
    group.by <- group.by %||% "ident"
    object[["ident"]] <- Idents(object = object)
    data <- object[[group.by]]
    for (group in group.by) {
      if (!is.factor(x = data[, group])) {
        data[, group] <- factor(x = data[, group])
      }
    }
  }
  else {
    if (interactive) {
      return(ISpatialFeaturePlot(object = object, feature = features[1], 
                                 image = images[1], slot = slot, alpha = alpha))
    }
    data <- FetchData(object = object, vars = features, slot = slot)
    features <- colnames(x = data)
    min.cutoff <- mapply(FUN = function(cutoff, feature) {
      return(ifelse(test = is.na(x = cutoff), yes = min(data[, 
                                                             feature]), no = cutoff))
    }, cutoff = min.cutoff, feature = features)
    max.cutoff <- mapply(FUN = function(cutoff, feature) {
      return(ifelse(test = is.na(x = cutoff), yes = max(data[, 
                                                             feature]), no = cutoff))
    }, cutoff = max.cutoff, feature = features)
    check.lengths <- unique(x = vapply(X = list(features, 
                                                min.cutoff, max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
    if (length(x = check.lengths) != 1) {
      stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    data <- sapply(X = 1:ncol(x = data), FUN = function(index) {
      data.feature <- as.vector(x = data[, index])
      min.use <- SetQuantile(cutoff = min.cutoff[index], 
                             data.feature)
      max.use <- SetQuantile(cutoff = max.cutoff[index], 
                             data.feature)
      data.feature[data.feature < min.use] <- min.use
      data.feature[data.feature > max.use] <- max.use
      return(data.feature)
    })
    colnames(x = data) <- features
    rownames(x = data) <- Cells(x = object)
  }
  features <- colnames(x = data)
  colnames(x = data) <- features
  rownames(x = data) <- colnames(x = object)
  facet.highlight <- facet.highlight && (!is.null(x = cells.highlight) && 
                                           is.list(x = cells.highlight))
  if (do.hover) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning("'do.hover' requires only one image, using image ", 
              images, call. = FALSE, immediate. = TRUE)
    }
    if (length(x = features) > 1) {
      features <- features[1]
      type <- ifelse(test = is.null(x = group.by), yes = "feature", 
                     no = "grouping")
      warning("'do.hover' requires only one ", type, ", using ", 
              features, call. = FALSE, immediate. = TRUE)
    }
    if (facet.highlight) {
      warning("'do.hover' requires no faceting highlighted cells", 
              call. = FALSE, immediate. = TRUE)
      facet.highlight <- FALSE
    }
  }
  if (facet.highlight) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning("Faceting the highlight only works with a single image, using image ", 
              images, call. = FALSE, immediate. = TRUE)
    }
    ncols <- length(x = cells.highlight)
  }
  else {
    ncols <- length(x = images)
  }
  plots <- vector(mode = "list", length = length(x = features) * 
                    ncols)
  for (i in 1:ncols) {
    plot.idx <- i
    image.idx <- ifelse(test = facet.highlight, yes = 1, 
                        no = i)
    image.use <- object[[images[[image.idx]]]]
    coordinates <- GetTissueCoordinates(object = image.use)
    highlight.use <- if (facet.highlight) {
      cells.highlight[i]
    }
    else {
      cells.highlight
    }
    for (j in 1:length(x = features)) {
      cols.unset <- is.factor(x = data[, features[j]]) && 
        is.null(x = cols)
      if (cols.unset) {
        cols <- hue_pal()(n = length(x = levels(x = data[, 
                                                         features[j]])))
        names(x = cols) <- levels(x = data[, features[j]])
      }
      plot <- SingleSpatialPlot(data = cbind(coordinates, 
                                             data[rownames(x = coordinates), features[j], 
                                                  drop = FALSE]), image = image.use, image.alpha = image.alpha, 
                                col.by = features[j], cols = cols, alpha.by = if (is.null(x = group.by)) {
                                  features[j]
                                }
                                else {
                                  NULL
                                }, pt.alpha = if (!is.null(x = group.by)) {
                                  alpha[j]
                                }
                                else {
                                  NULL
                                }, geom = if (inherits(x = image.use, what = "STARmap")) {
                                  "poly"
                                }
                                else {
                                  "spatial"
                                }, cells.highlight = highlight.use, cols.highlight = cols.highlight, 
                                pt.size.factor = pt.size.factor, stroke = stroke, 
                                crop = crop)
      if (is.null(x = group.by)) {
        plot <- plot + scale_fill_gradientn(name = features[j], 
                                            colours = SpatialColors(n = 100)) + theme(legend.position = "top") + 
          scale_alpha(range = alpha) + guides(alpha = FALSE)
      }
      else if (label) {
        plot <- LabelClusters(plot = plot, id = ifelse(test = is.null(x = cells.highlight), 
                                                       yes = features[j], no = "highlight"), geom = if (inherits(x = image.use, 
                                                                                                                 what = "STARmap")) {
                                                         "GeomPolygon"
                                                       }
                              else {
                                "GeomSpatial"
                              }, repel = repel, size = label.size, color = label.color, 
                              box = label.box, position = "nearest")
      }
      if (j == 1 && length(x = images) > 1 && !facet.highlight) {
        plot <- plot + ggtitle(label = images[[image.idx]]) + 
          theme(plot.title = element_text(hjust = 0.5))
      }
      if (facet.highlight) {
        plot <- plot + ggtitle(label = names(x = cells.highlight)[i]) + 
          theme(plot.title = element_text(hjust = 0.5)) + 
          NoLegend()
      }
      plots[[plot.idx]] <- plot
      plot.idx <- plot.idx + ncols
      if (cols.unset) {
        cols <- NULL
      }
    }
  }
  if (length(x = images) > 1 && combine) {
    plots <- wrap_plots(plots = plots, ncol = length(x = images))
  }
  else if (length(x = images == 1) && combine) {
    plots <- wrap_plots(plots = plots, ncol = ncol)
  }
  return(plots)
}

SpatialDimPlot <- function (object, group.by = NULL, images = NULL, cols = NULL, 
                            crop = TRUE, cells.highlight = NULL, cols.highlight = c("#DE2D26", 
                                                                                    "grey50"), facet.highlight = FALSE, label = FALSE, label.size = 7, 
                            label.color = "white", repel = FALSE, ncol = NULL, combine = TRUE, 
                            pt.size.factor = 1.6, alpha = c(1, 1), image.alpha = 1, stroke = 0.25, 
                            label.box = TRUE, interactive = FALSE, information = NULL) 
{
  return(SpatialPlot(object = object, group.by = group.by, 
                     images = images, cols = cols, crop = crop, cells.highlight = cells.highlight, 
                     cols.highlight = cols.highlight, facet.highlight = facet.highlight, 
                     label = label, label.size = label.size, label.color = label.color, 
                     repel = repel, ncol = ncol, combine = combine, pt.size.factor = pt.size.factor, 
                     alpha = alpha, image.alpha = image.alpha, stroke = stroke, 
                     label.box = label.box, interactive = interactive, information = information))
}


SpatialFeaturePlot <- function (object, features, images = NULL, crop = TRUE, slot = "data", 
                                min.cutoff = NA, max.cutoff = NA, ncol = NULL, combine = TRUE, 
                                pt.size.factor = 1.6, alpha = c(1, 1), image.alpha = 1, stroke = 0.25, 
                                interactive = FALSE, information = NULL) 
{
  return(SpatialPlot(object = object, features = features, 
                     images = images, crop = crop, slot = slot, min.cutoff = min.cutoff, 
                     max.cutoff = max.cutoff, ncol = ncol, combine = combine, 
                     pt.size.factor = pt.size.factor, alpha = alpha, image.alpha = image.alpha, 
                     stroke = stroke, interactive = interactive, information = information))
}
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))

circos.dendrogram2 <- function (dend, facing = c("outside", "inside"), max_height = NULL, 
                                use_x_attr = FALSE, sector.index = get.current.sector.index(), 
                                track.index = get.current.track.index(), xl=NULL) 
{
  os = get.current.sector.index()
  ot = get.current.track.index()
  set.current.cell(sector.index, track.index)
  on.exit(set.current.cell(os, ot))
  facing = match.arg(facing)[1]
  if (is.null(max_height)) {
    max_height = attr(dend, "height")
  }
  is.leaf = function(object) {
    leaf = attr(object, "leaf")
    if (is.null(leaf)) {
      FALSE
    }
    else {
      leaf
    }
  }
  use_x_attr = use_x_attr
  lines_par = function(col = par("col"), lty = par("lty"), 
                       lwd = par("lwd"), ...) {
    return(list(col = col, lty = lty, lwd = lwd))
  }
  points_par = function(col = par("col"), pch = par("pch"), 
                        cex = par("cex"), ...) {
    return(list(col = col, pch = pch, cex = cex))
  }
  draw.d = function(dend, max_height, facing = "outside", max_width = 0) {
    leaf = attr(dend, "leaf")
    height = attr(dend, "height")
    midpoint = attr(dend, "midpoint")
    n = length(dend)
    xl = numeric(n)
    yl = numeric(n)
    for (i in seq_len(n)) {
      if (use_x_attr) {
        xl[i] = attr(dend[[i]], "x")
      }
      else {
        if (is.leaf(dend[[i]])) {
          xl[i] = x[as.character(attr(dend[[i]], "label"))]
        }
        else {
          xl[i] = attr(dend[[i]], "midpoint") + x[as.character(labels(dend[[i]]))[1]]
        }
      }
      yl[i] = attr(dend[[i]], "height")
    }
    edge_par_lt = vector("list", n)
    for (i in seq_len(n)) {
      edge_par_lt[[i]] = do.call("lines_par", as.list(attr(dend[[i]], 
                                                           "edgePar")))
    }
    node_par = attr(dend, "nodePar")
    if (!is.null(node_par)) 
      node_par = do.call("points_par", as.list(attr(dend, 
                                                    "nodePar")))
    if (facing == "outside") {
      if (n == 1) {
        circos.lines(c(xl[1], xl[1]), max_height - c(yl[1], 
                                                     height), col = edge_par_lt[[1]]$col, lty = edge_par_lt[[1]]$lty, 
                     lwd = edge_par_lt[[1]]$lwd, straight = TRUE)
      }
      else {
        circos.lines(c(xl[1], xl[1]), max_height - c(yl[1], 
                                                     height), col = edge_par_lt[[1]]$col, lty = edge_par_lt[[1]]$lty, 
                     lwd = edge_par_lt[[1]]$lwd, straight = TRUE)
        circos.lines(c(xl[1], (xl[1] + xl[2])/2), max_height - 
                       c(height, height), col = edge_par_lt[[1]]$col, 
                     lty = edge_par_lt[[1]]$lty, lwd = edge_par_lt[[1]]$lwd)
        if (n > 2) {
          for (i in seq(2, n - 1)) {
            circos.lines(c(xl[i], xl[i]), max_height - 
                           c(yl[i], height), col = edge_par_lt[[i]]$col, 
                         lty = edge_par_lt[[i]]$lty, lwd = edge_par_lt[[i]]$lwd, 
                         straight = TRUE)
            circos.lines(c((xl[i - 1] + xl[i])/2, (xl[i] + 
                                                     xl[i + 1])/2), max_height - c(height, height), 
                         col = edge_par_lt[[i]]$col, lty = edge_par_lt[[i]]$lty, 
                         lwd = edge_par_lt[[i]]$lwd)
          }
        }
        circos.lines(c(xl[n], xl[n]), max_height - c(yl[n], 
                                                     height), col = edge_par_lt[[n]]$col, lty = edge_par_lt[[n]]$lty, 
                     lwd = edge_par_lt[[n]]$lwd, straight = TRUE)
        circos.lines(c(xl[n], (xl[n] + xl[n - 1])/2), 
                     max_height - c(height, height), col = edge_par_lt[[n]]$col, 
                     lty = edge_par_lt[[n]]$lty, lwd = edge_par_lt[[n]]$lwd)
      }
      if (!is.null(node_par)) {
        circos.points(mean(xl)/2, max_height - height, 
                      col = node_par$col, pch = node_par$pch, cex = node_par$cex)
      }
    }
    else if (facing == "inside") {
      if (n == 1) {
        circos.lines(c(xl[1], xl[1]), c(yl[1], height), 
                     col = edge_par_lt[[1]]$col, lty = edge_par_lt[[1]]$lty, 
                     lwd = edge_par_lt[[1]]$lwd, straight = TRUE)
      }
      else {
        circos.lines(c(xl[1], xl[1]), c(yl[1], height), 
                     col = edge_par_lt[[1]]$col, lty = edge_par_lt[[1]]$lty, 
                     lwd = edge_par_lt[[1]]$lwd, straight = TRUE)
        circos.lines(c(xl[1], (xl[1] + xl[2])/2), c(height, 
                                                    height), col = edge_par_lt[[1]]$col, lty = edge_par_lt[[1]]$lty, 
                     lwd = edge_par_lt[[1]]$lwd)
        if (n > 2) {
          for (i in seq(2, n - 1)) {
            circos.lines(c(xl[i], xl[i]), c(yl[i], height), 
                         col = edge_par_lt[[i]]$col, lty = edge_par_lt[[i]]$lty, 
                         lwd = edge_par_lt[[i]]$lwd, straight = TRUE)
            circos.lines(c((xl[i - 1] + xl[i])/2, (xl[i] + 
                                                     xl[i + 1])/2), c(height, height), col = edge_par_lt[[i]]$col, 
                         lty = edge_par_lt[[i]]$lty, lwd = edge_par_lt[[i]]$lwd)
          }
        }
        circos.lines(c(xl[n], xl[n]), c(yl[n], height), 
                     col = edge_par_lt[[n]]$col, lty = edge_par_lt[[n]]$lty, 
                     lwd = edge_par_lt[[n]]$lwd, straight = TRUE)
        circos.lines(c(xl[n], (xl[n] + xl[n - 1])/2), 
                     c(height, height), col = edge_par_lt[[n]]$col, 
                     lty = edge_par_lt[[n]]$lty, lwd = edge_par_lt[[n]]$lwd)
      }
      if (!is.null(node_par)) {
        circos.points(mean(xl)/2, height, col = node_par$col, 
                      pch = node_par$pch, cex = node_par$cex)
      }
    }
    for (i in seq_len(n)) {
      if (is.leaf(dend[[i]])) {
        node_par = attr(dend[[i]], "nodePar")
        if (!is.null(node_par)) 
          node_par = do.call("points_par", as.list(attr(dend[[i]], 
                                                        "nodePar")))
        if (facing == "outside") {
          if (!is.null(node_par)) {
            circos.points(xl[i], max_height, col = node_par$col, 
                          pch = node_par$pch, cex = node_par$cex)
          }
        }
        else if (facing == "inside") {
          if (!is.null(node_par)) {
            circos.points(xl[i], 0, col = node_par$col, 
                          pch = node_par$pch, cex = node_par$cex)
          }
        }
      }
      else {
        draw.d(dend[[i]], max_height, facing, max_width)
      }
    }
  }
  labels = as.character(labels(dend))
  x = xl
  names(x) = labels
  n = length(labels)
  if (!is.leaf(dend)) 
    draw.d(dend, max_height, facing, max_width = n)
}
