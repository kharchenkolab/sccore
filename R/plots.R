embeddingGroupPlot <- function(plot.df, groups, geom_point_w, min.cluster.size, mark.groups, font.size, legend.title, shuffle.colors, palette, ...) {
  groups <- as.factor(groups)

  plot.df$Group <- factor(NA, levels=levels(groups))
  arr.ids <- match(names(groups), plot.df$CellName)
  plot.df$Group[arr.ids[!is.na(arr.ids)]] <- groups[!is.na(arr.ids)]

  big.clusts <- (plot.df %>% subset(!is.na(Group)) %>% dplyr::group_by(Group) %>% dplyr::summarise(Size=n()) %>%
                   dplyr::filter(Size >= min.cluster.size))$Group %>% as.vector()

  plot.df$Group[!(plot.df$Group %in% big.clusts)] <- NA
  na.plot.df <- plot.df %>% dplyr::filter(is.na(Group))
  plot.df <- plot.df %>% dplyr::filter(!is.na(Group))

  gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x=x, y=y)) +
    geom_point_w(ggplot2::aes(col=Group))

  if (mark.groups) {
    labels.data <- plot.df %>% dplyr::group_by(Group) %>%
      dplyr::summarise(x=median(x), y=median(y), Size=n())

    if (length(font.size) == 1) {
      font.size <- c(font.size, font.size)
    }

    gg <- gg + ggrepel::geom_label_repel(
      data=labels.data, ggplot2::aes(label=Group, size=Size), color='black',
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

  color.vals <- palette(length(levels(groups)))
  if (shuffle.colors) {
    color.vals <- sample(color.vals)
  }
  gg <- gg + ggplot2::scale_color_manual(name=legend.title, values=color.vals, labels=levels(groups), drop=F) +
    ggplot2::guides(color=ggplot2::guide_legend(override.aes=list(alpha=1.0)))

  return(list(gg=gg, na.plot.df=na.plot.df))
}

embeddingColorsPlot <- function(plot.df, colors, groups, geom_point_w, gradient.range.quantile, color.range, legend.title, palette) {
  if(is.numeric(colors) && gradient.range.quantile < 1) {
    x <- colors;
    zlim <- as.numeric(quantile(x, p=c(1 - gradient.range.quantile, gradient.range.quantile), na.rm=TRUE))
    if(diff(zlim)==0) {
      zlim <- as.numeric(range(x))
    }
    x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
    colors <- x;
  }

  plot.df <- plot.df %>% dplyr::mutate(Color=colors[CellName])
  if(!is.null(groups)) {
    plot.df$Color[!plot.df$CellName %in% names(groups)] <- NA
  }
  na.plot.df <- plot.df %>% dplyr::filter(is.na(Color))
  plot.df <- plot.df %>% dplyr::filter(!is.na(Color))

  gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x=x, y=y))
  if(is.character(colors)) {
    gg <- gg + geom_point_w(color=plot.df$Color)
  } else {
    gg <- gg + geom_point_w(ggplot2::aes(col=Color))

    if (length(color.range) == 1) {
      if (color.range == "all") {
        color.range <- range(colors)
      } else if (color.range == "data") {
        color.range <- NULL
      } else {
        stop("Unknown color.range: ", color.range)
      }
    }

    if (!is.null(palette)) {
      gg <- gg + ggplot2::scale_colour_gradientn(colors=palette(100), limits=color.range)
    } else {
      c.range <- if (is.null(color.range)) range(colors, na.rm=T) else color.range
      if (prod(c.range) < 0) {
        gg <- gg + ggplot2::scale_color_gradient2(low="#0000ff",mid="#d8d0d0", high="#ff0000", limits=color.range)
      } else {
        gg <- gg + ggplot2::scale_color_gradient(low="#d8d0d0", high="#ff0000", limits=color.range)
      }
    }
  }

  if (!is.null(legend.title)) {
    gg <- gg + ggplot2::guides(color=ggplot2::guide_colorbar(title=legend.title))
  }

  return(list(gg=gg, na.plot.df=na.plot.df))
}

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
#' @param groups vector of cluster labels, names contain cell names
#' @param colors vector of numbers, which must be shouwn with point colors, names contain cell names. This argument is ignored if groups are provided.
#' @param subgroups subset of 'groups', selecting the cells for plot. Ignored if 'groups' is NULL
#' @param plot.na plot points, for which groups / colors are missed (TRUE / FALSE)
#' @param min.cluster.size labels for all groups with number of cells fewer than this parameter are considered as missed. This argument is ignored if groups aren't provided
#' @param mark.groups plot cluster labels above points
#' @param show.legend show legend
#' @param alpha opacity level [0; 1]
#' @param size point size
#' @param title plot title
#' @param plot.theme theme for the plot
#' @param palette function, which accepts number of colors and return list of colors (i.e. see colorRampPalette)
#' @param color.range controls range, in which colors are estimated. Pass "all" to estimate range based on all values of "colors", pass "data" to estimate it only based on colors, presented in the embedding. Alternatively you can pass vector of length 2 with (min, max) values.
#' @param font.size font size for cluster labels. It can either be single number for constant font size or pair (min, max) for font size depending on cluster size
#' @param show.ticks show ticks and tick labels
#' @param legend.position vector with (x, y) positions of the legend
#' @param legend.title legend title
#' @param gradient.range.quantile Winsorization quantile for the numeric colors and gene gradient
#' @param raster should layer with the points be rasterized (TRUE/ FALSE)? Setting of this argument to TRUE is useful when you need to export a plot with large number of points
#' @param raster.width width of the plot in inches. Ignored if raster == FALSE.
#' @param raster.height height of the plot in inches. Ignored if raster == FALSE.
#' @param raster.dpi dpi of the rasterized plot. Ignored if raster == FALSE.
#' @param shuffle.colors shuffle colors
#' @param keep.limits Keep axis limits from original plot, useful when plotting subgroups, only meaningful it plot.na=F
#' @return ggplot2 object
#' @export
embeddingPlot <- function(embedding, groups=NULL, colors=NULL, subgroups=NULL, plot.na=is.null(subgroups), min.cluster.size=0, mark.groups=TRUE,
                          show.legend=FALSE, alpha=0.4, size=0.8, title=NULL, plot.theme=NULL, palette=NULL, color.range="all",
                          font.size=c(3, 7), show.ticks=FALSE, show.labels=FALSE, legend.position=NULL, legend.title=NULL,
                          gradient.range.quantile=1, raster=FALSE, raster.width=NULL, raster.height=NULL, raster.dpi=300,
                          shuffle.colors=FALSE, keep.limits=!is.null(subgroups),
                          ...) {
  plot.df <- tibble::rownames_to_column(as.data.frame(embedding), "CellName")
  colnames(plot.df)[2:3] <- c("x", "y")

  if (raster && requireNamespace("ggrastr", quietly = TRUE)) {
    if (packageVersion("ggrastr") <= "0.1.6") {
      geom_point_w0 <- function(...)
        ggrastr::geom_point_rast(..., width=raster.width, height=raster.height, dpi=raster.dpi)
    } else {
      geom_point_w0 <- function(...)
        ggrastr::geom_point_rast(..., raster.width=raster.width, raster.height=raster.height, raster.dpi=raster.dpi)
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
    plot.info <- embeddingGroupPlot(plot.df, groups, geom_point_w, min.cluster.size, mark.groups,
                                    font.size, legend.title, shuffle.colors, palette, ...)
  } else if (!is.null(colors)) {
    plot.info <- embeddingColorsPlot(plot.df, colors, groups, geom_point_w, gradient.range.quantile,
                                     color.range, legend.title, palette)
  } else {
    plot.info <- list(gg=ggplot2::ggplot(plot.df, ggplot2::aes(x=x, y=y)) +
                        geom_point_w(alpha=alpha, size=size))
  }

  gg <- plot.info$gg
  if (plot.na && !is.null(plot.info$na.plot.df)) {
    gg <- gg + geom_point_w(data=plot.info$na.plot.df, color='black', shape=4)
  }

  if(keep.limits) {
    gg <- gg + ggplot2::lims(x=range(embedding[,1]), y=range(embedding[,2]))
  }

  gg <- styleEmbeddingPlot(gg, plot.theme=plot.theme, title=title, legend.position=legend.position,
                           show.legend=show.legend, show.ticks=show.ticks, show.labels=show.labels)
  return(gg)
}
