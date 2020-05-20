#' @import ggplot2
NULL

#' Factor to Color
#'
#' @description a utility function to translate factor into colors
#' @export
fac2col <- function(x, s=1, v=1, shuffle=FALSE, min.group.size=1,
                      return.details=FALSE, unclassified.cell.color='gray50', level.colors=NULL) {
  nx <- names(x)
  x <- as.factor(x)

  if (min.group.size>1) {
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

#' Encodes logic of how to handle named-vector and functional palettes
fac2palette <- function(groups, palette, unclassified.cell.color='gray50') {
  groups <- as.factor(groups)

  if (class(palette)=='function') {
    return(palette(length(levels(groups))))
  }

  if (is.list(palette)) { 
    palette <- setNames(unlist(palette),names(palette)) 
  }
  if (is.vector(palette)) {
    if (any(levels(groups) %in% names(palette))) {
      cols <- setNames(palette[match(levels(groups), names(palette))], levels(groups));
      cols[is.na(cols)] <- unclassified.cell.color
      return(cols)
    } else {
      # just take first n?
      if(length(palette)<length(levels(groups))) stop("provided palette does not have enough colors to show ",length(levels(groups))," levels")
      return(setNames(palette[1:length(levels(groups))],levels(groups)))
    }
  }
}

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

  color.vals <- fac2palette(groups,palette);


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
        color.range <- range(na.omit(colors))
      } else if (color.range == "data") {
        color.range <- NULL
      } else if (color.range == "symmetric") {
        if(is.numeric(colors)) {
          if(all(sign(na.omit(colors))>=0)) {
            color.range <- range(na.omit(colors))
          } else {
            color.range <- c(-1,1)*max(abs(colors))
          }
        } else {
          color.range <- NULL
        }
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
                          show.legend=FALSE, alpha=0.4, size=0.8, title=NULL, plot.theme=NULL, palette=NULL, color.range="symmetric",
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

  if (!is.null(subgroups) && !is.null(groups)) {
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

  if (keep.limits) {
    gg <- gg + ggplot2::lims(x=range(embedding[,1]), y=range(embedding[,2]))
  }

  gg <- styleEmbeddingPlot(gg, plot.theme=plot.theme, title=title, legend.position=legend.position,
                           show.legend=show.legend, show.ticks=show.ticks, show.labels=show.labels)
  return(gg)
}

#' Dot plot for group markers
#'
#' @description Dot plot adapted from Seurat:::DotPlot, see man for description.
#' @param markers Vector of gene markers to plot
#' @param count.matrix Merged count matrix, e.g., through conos.obj$getJointCountMatrix()
#' @param cell.groups Named factor containing cell groups (clusters) and cell names
#' @param marker.colour Character or numeric vector
#' @param cluster.colour Character or numeric vector
#' @param xlab X axis title
#' @param ylab Y axis title
#' @param ... Additional input to sccore:::plapply, see man for description.
#' @return ggplot2 object
#' @export
dotPlot <- function (markers,
                     count.matrix,
                     cell.groups,
                     verbose=TRUE,
                     n.cores=1,
                     marker.colour="black",
                     cluster.colour="black",
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
                     xlab = "Marker",
                     ylab = "Cluster",
                     ...) {
  
  scale.func <- switch(scale.by, 'size' = scale_size, 'radius' = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  if (!is.character(markers)) {
    stop("'markers' must be a character vector.")
  }

  missing.markers <- setdiff(markers, colnames(count.matrix))
  if (length(missing.markers)>0) {
    cat("Not all markers are in 'count.matrix'. The following are missing:\n",paste(missing.markers, collapse=" "),"\n")
    stop("Please update 'markers'.")
  }

  marker.table <- table(markers)
  if (sum(marker.table>1)!=0) {
    cat("The following genes are present more than once in 'markers':\n", paste(names(marker.table[marker.table>1]), collapse = " "), "\nThese genes will only be plotted at first instace. Consider revising.\n")
  }
  if (verbose) {
    cat("Extracting gene expression...\n")
  }
  # From CellAnnotatoR:::plotExpressionViolinMap, should be exchanged with generic function
  p.df <- plapply(markers, function(g) data.frame(Expr = count.matrix[names(cell.groups), g], Type = cell.groups, Gene = g), n.cores=n.cores, progress=verbose, ...) %>% Reduce(rbind, .)
  if (is.logical(gene.order) && gene.order) {
    gene.order <- unique(markers)
  } else {
    gene.order <- NULL
  }

  if (!is.null(gene.order)) {
    p.df %<>% dplyr::mutate(Gene = factor(as.character(Gene),
                                          levels = gene.order))
  }

  # Adapted from Seurat:::DotPlot
  if (verbose) { 
    cat("Calculating expression distributions...\n")
  }
  data.plot <- levels(cell.groups) %>% plapply(function(t) {
    markers %>% lapply(function(g) {
      df <- p.df %>% filter(Type==t, Gene==g)
      pct.exp <- sum(df$Expr>0)/dim(df)[1]*100
      avg.exp <- mean(df$Expr[df$Expr>0])
      res <- data.frame(gene=g,
                        pct.exp=pct.exp,
                        avg.exp=avg.exp)
      return(res)
    }) %>% Reduce(rbind, .)
  }, n.cores=n.cores, progress=verbose, ...) %>%
    setNames(levels(cell.groups)) %>%
    dplyr::bind_rows(., .id="cluster")

  data.plot$cluster %<>% factor(., levels=rev(unique(.)))

  data.plot %<>% dplyr::arrange(gene)

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
