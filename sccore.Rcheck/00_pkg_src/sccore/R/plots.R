#' @import ggplot2
#' @importFrom graphics par
#' @importFrom grDevices adjustcolor rainbow colorRampPalette
#' @importFrom magrittr %>% %<>% %$%
#' @importFrom rlang .data
#' @import scales
NULL

## for magrittr and dplyr functions below
if(getRversion() >= "2.15.1"){
  utils::globalVariables(c(".", "n", "Size", "Group", "x", "y"))
}


#' Utility function to translate a factor into colors
#'
#' @param x input factor
#' @param s numeric The "saturation" to be used to complete the HSV color descriptions (default=1) See ?rainbow in Palettes, grDevices
#' @param v numeric The "value" to be used to complete the HSV color descriptions (default=1) See ?rainbow in Palettes, grDevices
#' @param shuffle boolean If TRUE, shuffles columns with shuffle(columns) (default=FALSE)
#' @param min.group.size integer Exclude groups of size less than the min.group.size (default=1)
#' @param return.details boolean If TRUE, returns a list list(colors=y, palette=col). Otherwise, just returns the factor (default=FALSE)
#' @param unclassified.cell.color Color for unclassified cells (default='gray50')
#' @param level.colors (default=NULL)
#' @return vector or list of colors
#' @examples
#' genes = factor(c("BRAF", "NPC1", "PAX3", "BRCA2", "FMR1"))
#' fac2col(genes)
#'
#' @export
fac2col <- function(x, s=1, v=1, shuffle=FALSE, min.group.size=1,
                      return.details=FALSE, unclassified.cell.color='gray50', level.colors=NULL) {
  nx <- names(x)
  x <- as.factor(x)

  if (min.group.size > 1) {
    x <- factor(x, exclude=levels(x)[unlist(tapply(rep(1,length(x)), x, length)) < min.group.size])
    x <- droplevels(x)
  }

  if (is.null(level.colors)) {
    col <- rainbow(length(levels(x)), s=s, v=v)
  } else {
    col <- level.colors[1:length(levels(x))]
  }

  names(col) <- levels(x)

  if (shuffle){
    col <- sample(col)
  }

  y <- col[as.integer(x)]
  names(y) <- names(x)
  y[is.na(y)] <- unclassified.cell.color
  names(y) <- nx

  if (return.details) {
    return(list(colors=y, palette=col))
  } else {
    return(y)
  }
}

#' Utility function to translate values into colors.
#'
#' @param x input values
#' @param gradientPalette gradient palette (default=NULL). If NULL, use colorRampPalette(c('gray90','red'), space = "Lab")(1024) if the values are non-negative; otherwise colorRampPalette(c("blue", "grey90", "red"), space = "Lab")(1024) is used
#' @param zlim a two-value vector specifying limits of the values that should correspond to the extremes of the color gradient
#' @param gradient.range.quantile extreme quantiles of values that should be trimmed prior to color mapping (default=0.95)
#' @examples
#' colors <- val2col( rnorm(10) )
#' 
#' @export
val2col <- function(x, gradientPalette=NULL, zlim=NULL, gradient.range.quantile=0.95) {
  nx <- names(x);
  if (all(sign(x)>=0)) {
    if(is.null(gradientPalette)) {
      gradientPalette <- colorRampPalette(c('gray90','red'), space = "Lab")(1024)
    }
    if (is.null(zlim)) {
      zlim <- as.numeric(quantile(na.omit(x),p=c(1-gradient.range.quantile,gradient.range.quantile)))
      if(diff(zlim)==0) {
        zlim <- as.numeric(range(na.omit(x)))
      }
    }
    x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
    x <- (x-zlim[1])/(zlim[2]-zlim[1])

  } else {
    if(is.null(gradientPalette)) {
      gradientPalette <- colorRampPalette(c("blue", "grey90", "red"), space = "Lab")(1024)
    }
    if(is.null(zlim)) {
      zlim <- c(-1,1)*as.numeric(quantile(na.omit(abs(x)),p=gradient.range.quantile))
      if(diff(zlim)==0) {
        zlim <- c(-1,1)*as.numeric(na.omit(max(abs(x))))
      }
    }
    x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
    x <- (x-zlim[1])/(zlim[2]-zlim[1])

  }

  col <- gradientPalette[x*(length(gradientPalette)-1)+1]
  names(col) <- nx
  return(col)
}


#' Encodes logic of how to handle named-vector and functional palettes. Used primarily within embeddingGroupPlot()
#'
#' @param groups vector of cluster labels, names contain cell names
#' @param palette function, which accepts number of colors and return list of colors (i.e. see 'colorRampPalette')
#' @param unclassified.cell.color Color for unclassified cells (default='gray50')
#' @return vector or palette
fac2palette <- function(groups, palette, unclassified.cell.color='gray50') {
  groups <- as.factor(groups)

  if (class(palette)=='function') {
    return(palette(length(levels(groups))))
  }

  if (is.list(palette)) {
    palette <- stats::setNames(unlist(palette),names(palette))
  }
  if (is.vector(palette)) {
    if (any(levels(groups) %in% names(palette))) {
      cols <- stats::setNames(palette[match(levels(groups), names(palette))], levels(groups));
      cols[is.na(cols)] <- unclassified.cell.color
      return(cols)
    } else {
      # just take first n?
      if(length(palette)<length(levels(groups))) stop("provided palette does not have enough colors to show ",length(levels(groups))," levels")
      return(stats::setNames(palette[1:length(levels(groups))],levels(groups)))
    }
  }
}

#' Helper function to return a ggplot color gradient for a numeric vector
#' ggplot(aes(color=x, ...), ...) + val2ggcol(x)
#'
#' @param values values by which the color gradient is determined
#' @param gradient.range.quantile numeric Trimming quantile (default=1). Either a single number or two numbers - for lower and upper quantile.
#' @param color.range either a vector of two values explicitly specifying the values corresponding to the start/end of the gradient, or string "symmetric" or "all" (default="symmetric"). "symmetric": range will fit data, but will be symmetrized around zeros, "all": gradient will match the span of the range of the data (after gradient.range.quantile)
#' @param palette an optional palette fucntion (default=NULL). The default becomes blue-gray90-red; if the values do not straddle 0, then truncated gradients (blue-gray90 or gray90-red) will be used
#' @param midpoint optional midpoint (default=NULL). Set for the center of the resulting range by default
#' @param oob function to determine what to do with the values outside of the range (default =scales::squish). Refer to 'oob' parameter in ggplot
#' @param return.fill boolean Whether to return fill gradients instead of color (default=FALSE)
#' @param ... additional arguments are passed to ggplot2::scale_color_gradient* functions, i.e. scale_color_gradient(), scale_color_gradient2(), scale_color_gradientn()
#' @return ggplot2::scale_colour_gradient object
val2ggcol <- function(values, gradient.range.quantile=1, color.range='symmetric', palette=NULL, midpoint=NULL, oob=scales::squish, return.fill=FALSE, ...) {
  if(length(gradient.range.quantile)>1) { # min/max quantile is given
    zlim <- as.numeric(quantile(values, p=gradient.range.quantile, na.rm=TRUE))
  } else if(gradient.range.quantile<1) { # single value spec
    zlim <- sort(as.numeric(quantile(values, p=c(1 - gradient.range.quantile, gradient.range.quantile), na.rm=TRUE)))
  } else {
    zlim <- range(stats::na.omit(values))
  }

  ## Symmetrize the range for vectors that span 0.
  ## Vectors that are squarely in the positive or negative territory are not symmetrized.
  if (length(color.range) == 1) {
    if (!(color.range %in% c('symmetric', 'all'))) {
      stop("Can't parse color.range: ", color.range)
    }

    if((color.range == 'symmetric') && (prod(zlim) < 0)) {
      zlim <- c(-1, 1)*max(abs(zlim))
    }
  } else if (length(color.range) == 2) {
    zlim <- color.range
  } else {
    stop("Can't parse color.range: ", color.range)
  }

  if(is.null(midpoint)){
    midpoint <- sum(zlim)/2
  }

  # pick a palette and return
  if (is.null(palette)) {
    if (max(abs(zlim))==0) {
      ## if gene counts all 0, then simply plot all cells as "gray90"
      ggplot2::scale_color_gradient(low="gray90", high="gray90", limits=zlim, ...)
    } else if(zlim[2]<0) {
      if(return.fill) {
        ggplot2::scale_fill_gradient(low="blue", high="gray90", limits=zlim, oob=oob, ...)
      } else {
        ggplot2::scale_color_gradient(low="blue", high="gray90", limits=zlim, oob=oob, ...)
      }
    } else if(zlim[1]>=0) {
      if(return.fill) {
        ggplot2::scale_fill_gradient(low="gray90", high="red", limits=zlim, oob=oob,  ...)
      } else {
        ggplot2::scale_color_gradient(low="gray90", high="red", limits=zlim, oob=oob,  ...)
      }
    } else {
      if(return.fill) {
        ggplot2::scale_fill_gradient2(low="blue",mid="gray90", high="red", limits=zlim, oob=oob, midpoint=midpoint, ...)
      } else {
        ggplot2::scale_color_gradient2(low="blue",mid="gray90", high="red", limits=zlim, oob=oob, midpoint=midpoint, ...)
      }
    }
  } else {
    if(return.fill) {
      ggplot2::scale_fill_gradientn(colors=palette(100), limits=zlim, oob=oob, ...)
    } else {
      ggplot2::scale_colour_gradientn(colors=palette(100), limits=zlim, oob=oob, ...)
    }
  }
}

#' Plotting function for cluster labels, names contain cell names. Used primarily in embeddingPlot().
#'
#' @inheritParams embeddingPlot
#' @param plot.df data.frame for plotting. In embeddingPlot(), this is a tibble from tibble::rownames_to_column().
#' @param geom_point_w function to work with geom_point layer from ggplot2 (default=ggplot2::geom_point)
#' @param ... Additional arguments passed to ggplot2::geom_label_repel()
#' @return ggplot2 object
embeddingGroupPlot <- function(plot.df, groups, geom_point_w, min.cluster.size, mark.groups, font.size, legend.title, shuffle.colors, palette, plot.na, ...) {

  groups <- as.factor(groups)

  plot.df$Group <- factor(NA, levels=levels(groups))
  arr.ids <- match(names(groups), plot.df$CellName)
  plot.df$Group[arr.ids[!is.na(arr.ids)]] <- groups[!is.na(arr.ids)]

  big.clusts <- (plot.df %>% subset(!is.na(Group)) %>% dplyr::group_by(Group) %>% dplyr::summarise(Size=dplyr::n()) %>%
                   dplyr::filter(Size >= min.cluster.size))$Group %>% as.vector()

  plot.df$Group[!(plot.df$Group %in% big.clusts)] <- NA
  na.plot.df <- plot.df %>% dplyr::filter(is.na(Group))
  plot.df <- plot.df %>% dplyr::filter(!is.na(Group))

  gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x=x, y=y))

  ## If plot.na passed a numeric value below 0, the NA symbols are plotted below the cells.
  ## Otherwise they’re plotted above the cells.
  if (plot.na & (plot.na < 0)) {
    gg <- gg + geom_point_w(data=na.plot.df, color='black', shape=4)
  }

  gg <- gg + geom_point_w(ggplot2::aes(col=.data$Group))

  if (plot.na & (plot.na >= 0)) {
    gg <- gg + geom_point_w(data=na.plot.df, color='black', shape=4)
  }

  if (mark.groups) {
    labels.data <- plot.df %>% dplyr::group_by(Group) %>%
      dplyr::summarise(x=median(x), y=median(y), Size=dplyr::n())

    if (length(font.size) == 1) {
      font.size <- c(font.size, font.size)
    }

    gg <- gg + ggrepel::geom_label_repel(
      data=labels.data, ggplot2::aes(label=.data$Group, size=.data$Size), color='black',
      fill=ggplot2::alpha('white', 0.7), label.size = NA,
      label.padding=ggplot2::unit(1, "pt"), seed=42, ...) +
      ggplot2::scale_size_continuous(range=font.size, trans='identity', guide='none')
  }

  if (is.null(legend.title)) {
    legend.title <- "Group"
  }

  if(is.null(palette)) {
    palette <- rainbow
  }

  color.vals <- fac2palette(groups,palette);


  if (shuffle.colors) {
    color.vals <- sample(color.vals)
  }
  gg <- gg + ggplot2::scale_color_manual(name=legend.title, values=color.vals, labels=levels(groups), drop=FALSE) +
    ggplot2::guides(color=ggplot2::guide_legend(override.aes=list(alpha=1.0)))

  return(gg)
}

#' Set colors for embedding plot. Used primarily in embeddingPlot().
#'
#' @inheritParams embeddingPlot
#' @param plot.df data.frame for plotting. In embeddingPlot(), this is a tibble from tibble::rownames_to_column().
#' @param geom_point_w function to work with geom_point layer from ggplot2 (default=ggplot2::geom_point)
#' @return ggplot2 object
embeddingColorsPlot <- function(plot.df, colors, groups=NULL, geom_point_w=ggplot2::geom_point, gradient.range.quantile=1, color.range="symmetric", legend.title=NULL, palette=NULL, plot.na=TRUE) {
  plot.df <- plot.df %>% dplyr::mutate(Color=colors[.data$CellName])
  if(!is.null(groups)) {
    plot.df$Color[!plot.df$CellName %in% names(groups)] <- NA
  }
  na.plot.df <- plot.df %>% dplyr::filter(is.na(.data$Color))
  plot.df <- plot.df %>% dplyr::filter(!is.na(.data$Color))

  gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x=x, y=y))

  if (plot.na & (plot.na < 0)) {
    gg <- gg + geom_point_w(data=na.plot.df, color='black', shape=4)
  }

  if(is.character(colors)) {
    gg <- gg + geom_point_w(color=plot.df$Color)
  } else {
    gg <- gg + geom_point_w(ggplot2::aes(col=.data$Color)) + val2ggcol(plot.df$Color, gradient.range.quantile=gradient.range.quantile, palette=palette, color.range=color.range)
  }

  if (plot.na & (plot.na >= 0)) {
    gg <- gg + geom_point_w(data=na.plot.df, color='black', shape=4)
  }

  if (!is.null(legend.title)) {
    gg <- gg + ggplot2::guides(color=ggplot2::guide_colorbar(title=legend.title))
  }

  return(gg)
}

#' Set plot.theme, legend, ticks for embedding plot. Used primarily in embeddingPlot().
#'
#' @param gg ggplot2 object to plot
#' @param plot.theme theme for the plot (default=NULL)
#' @param title plot title (default=NULL)
#' @param legend.position vector with (x, y) positions of the legend (default=NULL)
#' @param show.legend show legend (default=TRUE)
#' @param show.ticks show ticks and tick labels (default=TRUE)
#' @param show.labels show labels (default=TRUE)
#' @param relabel.axis boolean If TRUE, relabel axes with ggplot2::labs(x='Component 1', y='Component 2') (default=TRUE)
#' @return ggplot2 object
styleEmbeddingPlot <- function(gg, plot.theme=NULL, title=NULL, legend.position=NULL, show.legend=TRUE, show.ticks=TRUE, show.labels=TRUE, relabel.axis=TRUE) {
  if (relabel.axis) {
    gg <- gg + ggplot2::labs(x='Component 1', y='Component 2')
  }

  if (!is.null(plot.theme)) {
    gg <- gg + plot.theme
  }

  if (!is.null(title)) {
    gg <- gg + ggplot2::ggtitle(title)
  }

  if (!is.null(legend.position)) {
    gg <- gg + ggplot2::theme(legend.position=legend.position,
                              legend.justification=legend.position)
  }

  if (!show.legend) {
    gg <- gg + ggplot2::theme(legend.position="none")
  }

  if (!show.ticks) {
    gg <- gg + ggplot2::theme(axis.ticks=ggplot2::element_blank(),
                              axis.text=ggplot2::element_blank())
  }

  if (!show.labels) {
    gg <- gg + ggplot2::theme(axis.title=ggplot2::element_blank())
  }

  return(gg)
}

#' Plot embedding with provided labels / colors using ggplot2
#'
#' @inheritDotParams ggrepel::geom_label_repel
#' @param embedding two-column matrix with x and y coordinates of the embedding, rownames contain cell names and are used to match coordinates with groups or colors
#' @param groups vector of cluster labels, names contain cell names (default=NULL)
#' @param colors vector of numbers, which must be shown with point colors, names contain cell names (default=NULL). This argument is ignored if groups are provided.
#' @param subgroups subset of 'groups', selecting the cells for plot (default=NULL). Ignored if 'groups' is NULL
#' @param plot.na boolean/numeric Whether to plot points, for which groups / colors are missed (default=is.null(subgroups), i.e. FALSE). If plot.na passed a numeric value below 0, the NA symbols are plotted below the cells. Otherwise if values >=0, they’re plotted above the cells. Note that this argument is FALSE if 'subgroups' is NULL
#' @param min.cluster.size labels for all groups with number of cells fewer than this parameter are considered as missed (default=0). This argument is ignored if groups aren't provided
#' @param mark.groups plot cluster labels above points (default=TRUE)
#' @param show.legend show legend (default=FALSE)
#' @param alpha opacity level [0, 1] (default=0.4)
#' @param size point size (default=0.8)
#' @param title plot title (default=NULL)
#' @param plot.theme theme for the plot (default=NULL)
#' @param palette function, which accepts number of colors and return list of colors (i.e. see 'colorRampPalette') (default=NULL)
#' @param color.range controls range, in which colors are estimated (default="symmetric"). Pass "all" to estimate range based on all values of "colors", pass "data" to estimate it only based on colors, presented in the embedding. Alternatively you can pass vector of length 2 with (min, max) values.
#' @param font.size font size for cluster labels (default=c(3, 7)). It can either be single number for constant font size or pair (min, max) for font size depending on cluster size
#' @param show.ticks show ticks and tick labels (default=FALSE)
#' @param show.labels show labels (default=FALSE)
#' @param legend.position vector with (x, y) positions of the legend (default=NULL)
#' @param legend.title legend title (default=NULL)
#' @param gradient.range.quantile Winsorization quantile for the numeric colors and gene gradient (default=1)
#' @param raster boolean whether layer with the points be rasterized (default=FALSE). Setting of this argument to TRUE is useful when you need to export a plot with large number of points
#' @param raster.dpi dpi of the rasterized plot. (default=300). Ignored if raster == FALSE.
#' @param shuffle.colors shuffle colors (default=FALSE)
#' @param keep.limits Keep axis limits from original plot (default=!is.null(subgroups)). Useful when plotting subgroups, only meaningful it plot.na=FALSE
#' @return ggplot2 object
#' @examples
#' library(sccore)
#' embeddingPlot(umapEmbedding, show.ticks=TRUE, show.labels=TRUE, title="UMAP embedding")
#'
#' @export
embeddingPlot <- function(embedding, groups=NULL, colors=NULL, subgroups=NULL, plot.na=is.null(subgroups), min.cluster.size=0, mark.groups=TRUE,
                          show.legend=FALSE, alpha=0.4, size=0.8, title=NULL, plot.theme=NULL, palette=NULL, color.range="symmetric",
                          font.size=c(3, 7), show.ticks=FALSE, show.labels=FALSE, legend.position=NULL, legend.title=NULL,
                          gradient.range.quantile=1, raster=FALSE, raster.dpi=300, shuffle.colors=FALSE, keep.limits=!is.null(subgroups), ...) {
  plot.df <- tibble::rownames_to_column(as.data.frame(embedding), "CellName")
  colnames(plot.df)[2:3] <- c("x", "y")

  if (raster && requireNamespace("ggrastr", quietly = TRUE)) {
    if (utils::packageVersion("ggrastr") <= "0.1.6") {
      geom_point_w0 <- function(...)
        ggrastr::geom_point_rast(..., dpi=raster.dpi)
    } else {
      geom_point_w0 <- function(...)
        ggrastr::geom_point_rast(..., raster.dpi=raster.dpi)
    }
  } else {
    if (raster) {
      warning("You have to install ggrastr package to be able to use 'raster' parameter")
    }
    geom_point_w0 <- ggplot2::geom_point
  }

  geom_point_w <- function(...) geom_point_w0(..., size=size, alpha=alpha)

  if(!is.null(subgroups) && !is.null(groups)) {
    groups %<>% .[. %in% subgroups]
    if(length(groups)==0) {
      stop("'groups' is empty after filtering by 'subgroups'.")
    }
  }

  if (!is.null(groups) && is.null(colors)) {
    gg <- embeddingGroupPlot(plot.df, groups, geom_point_w, min.cluster.size, mark.groups, font.size,
                             legend.title, shuffle.colors, palette, plot.na=plot.na, ...)
  } else if (!is.null(colors)) {
    gg <- embeddingColorsPlot(plot.df, colors, groups, geom_point_w, gradient.range.quantile,
                              color.range, legend.title, palette, plot.na=plot.na)
  } else {
    gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x=x, y=y)) + geom_point_w()
  }

  if(keep.limits) {
    gg <- gg + ggplot2::lims(x=range(embedding[,1]), y=range(embedding[,2]))
  }

  gg <- styleEmbeddingPlot(gg, plot.theme=plot.theme, title=title, legend.position=legend.position,
                           show.legend=show.legend, show.ticks=show.ticks, show.labels=show.labels)
  return(gg)
}

#' Dot plot adapted from Seurat:::DotPlot, see ?Seurat:::DotPlot for details
#'
#' @param markers Vector of gene markers to plot
#' @param count.matrix Merged count matrix, cells in rows and genes in columns
#' @param cell.groups Named factor containing cell groups (clusters) and cell names as names
#' @param marker.colour Character or numeric vector (default="black")
#' @param cluster.colour Character or numeric vector (default="black")
#' @param xlab string X-axis title (default="Marker")
#' @param ylab string Y-axis title (default="Cluster")
#' @param n.cores integer Number of cores (default=1)
#' @param text.angle numeric Angle of text displayed (default=45)
#' @param gene.order Either factor of genes passed to dplyr::mutate(levels=gene.order), or a boolean. (default=NULL) If TRUE, gene.order is set to the unique markers. If FALSE, gene.order is set to NULL. If NULL, the argument is ignored.
#' @param cols Colors to plot (default=c("blue", "red")). The name of a palette from 'RColorBrewer::brewer.pal.info', a pair of colors defining a gradient, or 3+ colors defining multiple gradients (if 'split.by' is set).
#' @param col.min numeric Minimum scaled average expression threshold (default=-2.5). Everything smaller will be set to this.
#' @param col.max numeric Maximum scaled average expression threshold (default=2.5). Everything larger will be set to this.
#' @param dot.min numeric The fraction of cells at which to draw the smallest dot (default=0). All cell groups with less than this expressing the given gene will have no dot drawn.
#' @param dot.scale numeric Scale the size of the points, similar to cex (default=6)
#' @param scale.by  string Scale the size of the points by 'size' or by 'radius' (default="radius")
#' @param scale.min numeric Set lower limit for scaling, use NA for default (default=NA)
#' @param scale.max numeric Set upper limit for scaling, use NA for default (default=NA)
#' @param verbose boolean Verbose output (default=TRUE)
#' @param ... Additional inputs passed to sccore::plapply(), see man for description.
#' @return ggplot2 object
#' @examples
#' library(dplyr)
#' ## Create merged count matrix
#' ## In this example, cms is a list of count matrices from, e.g., Cellranger count, 
#' ## where cells are in columns and genes in rows
#' ## cm <- sccore:::mergeCountMatrices(cms, transposed = FALSE) %>% Matrix::t()
#'
#' ## If coming from Conos, this can be extracted like so
#' ## cm <- conos.obj$getJointCountMatrix(raw = FALSE) # Either normalized or raw values can be used
#'
#' ## Here, we create a random sparse matrix
#' cm <- Matrix::rsparsematrix(30,3,0.5) %>% abs(.) %>% 
#'             `dimnames<-`(list(1:30,c("gene1","gene2","gene3")))
#'
#' ## Create marker vector
#' markers <- c("gene1","gene2","gene3")
#'
#' ## Additionally, color vectors can be included. 
#' ## These should have the same length as the input (markers, cell groups) 
#' ## Otherwise, they are recycled
#' col.markers <- c("black","black","red") # or c(1,1,2)
#' col.clusters <- c("black","red","black") # or c(1,2,1)
#'
#' ## Create annotation vector
#' annotation <- c(rep("cluster1",10),rep("cluster2",10),rep("cluster3",10)) %>% 
#'     factor() %>% setNames(1:30)
#'
#' ## Plot. Here, the expression colours range from gray (low expression) to purple (high expression)
#' sccore:::dotPlot(markers = markers, count.matrix = cm, cell.groups = annotation, 
#'     marker.colour = col.markers, cluster.colour = col.clusters, cols=c("gray","purple"))
#'
#' @export
dotPlot <- function (markers,
                     count.matrix,
                     cell.groups,
                     marker.colour="black",
                     cluster.colour="black",
                     xlab = "Marker",
                     ylab = "Cluster",
                     n.cores = 1,
                     text.angle = 45,
                     gene.order = NULL,
                     cols = c("blue", "red"),
                     col.min = -2.5,
                     col.max = 2.5,
                     dot.min = 0,
                     dot.scale = 6,
                     scale.by = "radius",
                     scale.min = NA,
                     scale.max = NA,
                     verbose=TRUE,
                     ...) {

  scale.func <- switch(scale.by, 'size' = scale_size, 'radius' = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  if (!is.character(markers)) {
    stop("'markers' must be a character vector.")
  }

  missing.markers <- setdiff(markers, colnames(count.matrix))
  if (length(missing.markers)>0) {
    message("Not all markers are in 'count.matrix'. The following are missing: ",paste(missing.markers, collapse=" "))
    stop("Please update 'markers'.")
  }

  marker.table <- table(markers)
  if (sum(marker.table>1)!=0) {
    message("The following genes are present more than once in 'markers': ", paste(names(marker.table[marker.table>1]), collapse = " "), " These genes will only be plotted at first instace. Consider revising. ")
  }
  if (verbose) {
    message("Extracting gene expression... ")
  }

  if (class(cell.groups) != "factor") {
    tryCatch({
      if(verbose){
        message("Treating 'cell.groups' as a factor.")
      }
      cell.groups %<>% as.factor()
    }, error=function(e) stop("Could not convert 'cell.groups' to a factor\n", e))
  }
  # From CellAnnotatoR:::plotExpressionViolinMap, should be exchanged with generic function
  p.df <- plapply(markers, function(g) data.frame(Expr = count.matrix[names(cell.groups), g], Type = cell.groups, Gene = g), n.cores=n.cores, progress=verbose, ...) %>% Reduce(rbind, .)
  if (is.logical(gene.order) && gene.order) {
    gene.order <- unique(markers)
  } else {
    gene.order <- NULL
  }

  if (!is.null(gene.order)) {
    p.df %<>% dplyr::mutate(Gene = factor(as.character(.data$Gene), levels = gene.order))
  }

  # Adapted from Seurat:::DotPlot
  if (verbose) {
    message("Calculating expression distributions... ")
  }
  data.plot <- levels(cell.groups) %>% plapply(function(t) {
    markers %>% lapply(function(g) {
      df <- p.df %>% dplyr::filter(.data$Type==t, .data$Gene==g)
      pct.exp <- sum(df$Expr>0)/dim(df)[1]*100
      avg.exp <- mean(df$Expr[df$Expr>0])
      res <- data.frame(gene=g,
                        pct.exp=pct.exp,
                        avg.exp=avg.exp)
      return(res)
    }) %>% Reduce(rbind, .)
  }, n.cores=n.cores, progress=verbose, ...) %>%
    stats::setNames(levels(cell.groups)) %>%
    dplyr::bind_rows(., .id="cluster")

  data.plot$cluster %<>% factor(., levels=rev(unique(.)))

  data.plot %<>% dplyr::arrange(.data$gene)

  data.plot$avg.exp.scaled <- data.plot$gene %>% unique %>% sapply(function(g) {
    data.plot %>% .[.$gene == g, 'avg.exp'] %>%
      scale %>%
      setMinMax(min = col.min, max = col.max)
  }) %>% unlist %>% as.numeric

  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA

  cluster.colour %<>% rev

  plot <- ggplot(data.plot, aes_string("gene", "cluster")) +
    geom_point(aes_string(size = "pct.exp", color = "avg.exp.scaled")) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.text.x = element_text(angle=text.angle, hjust = 1, colour=marker.colour),
          axis.text.y = element_text(colour=cluster.colour),
          panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(size = guide_legend(title = 'Percent expressed'), color = guide_colorbar(title = 'Average expression')) +
    labs(x = xlab, y = ylab) +
    scale_color_gradient(low = cols[1], high = cols[2])

  return(plot)
}
