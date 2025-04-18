% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\docType{methods}
\name{embeddingPlot}
\alias{embeddingPlot}
\alias{embeddingPlot,ANY-method}
\alias{embeddingPlot-methods,ANY}
\alias{embeddingPlot,Seurat-method}
\alias{embeddingPlot-methods,Seurat}
\title{embeddingPlot generic
Plot embedding with provided labels / colors using ggplot2}
\usage{
embeddingPlot(object, ...)

\S4method{embeddingPlot}{ANY}(
  object,
  groups = NULL,
  colors = NULL,
  subgroups = NULL,
  plot.na = is.null(subgroups),
  min.cluster.size = 0,
  mark.groups = TRUE,
  show.legend = FALSE,
  alpha = 0.4,
  size = 0.8,
  title = NULL,
  plot.theme = NULL,
  palette = NULL,
  color.range = "symmetric",
  font.size = c(3, 7),
  show.ticks = FALSE,
  show.labels = FALSE,
  legend.position = NULL,
  legend.title = NULL,
  gradient.range.quantile = 1,
  raster = FALSE,
  raster.dpi = 300,
  shuffle.colors = FALSE,
  keep.limits = !is.null(subgroups),
  ...
)

\S4method{embeddingPlot}{Seurat}(object, reduction = NULL, groups = NULL, colors = NULL, ...)
}
\arguments{
\item{object}{Seurat object}

\item{...}{
  Arguments passed on to \code{\link[ggrepel:geom_text_repel]{ggrepel::geom_label_repel}}
  \describe{
    \item{\code{mapping}}{Set of aesthetic mappings created by \code{\link[ggplot2]{aes}} or
\code{\link[ggplot2]{aes_}}. If specified and \code{inherit.aes = TRUE} (the
default), is combined with the default mapping at the top level of the
plot. You only need to supply \code{mapping} if there isn't a mapping
defined for the plot.}
    \item{\code{data}}{A data frame. If specified, overrides the default data frame
defined at the top level of the plot.}
    \item{\code{stat}}{The statistical transformation to use on the data for this
layer, as a string.}
    \item{\code{position}}{Position adjustment, either as a string, or the result of
a call to a position adjustment function.}
    \item{\code{parse}}{If TRUE, the labels will be parsed into expressions and
displayed as described in ?plotmath}
    \item{\code{box.padding}}{Amount of padding around bounding box, as unit or number.
Defaults to 0.25. (Default unit is lines, but other units can be specified
by passing \code{unit(x, "units")}).}
    \item{\code{label.padding}}{Amount of padding around label, as unit or number.
Defaults to 0.25. (Default unit is lines, but other units can be specified
by passing \code{unit(x, "units")}).}
    \item{\code{point.padding}}{Amount of padding around labeled point, as unit or
number. Defaults to 0. (Default unit is lines, but other units can be
specified by passing \code{unit(x, "units")}).}
    \item{\code{label.r}}{Radius of rounded corners, as unit or number. Defaults
to 0.15. (Default unit is lines, but other units can be specified by
passing \code{unit(x, "units")}).}
    \item{\code{label.size}}{Size of label border, in mm.}
    \item{\code{min.segment.length}}{Skip drawing segments shorter than this, as unit or
number. Defaults to 0.5. (Default unit is lines, but other units can be
specified by passing \code{unit(x, "units")}).}
    \item{\code{arrow}}{specification for arrow heads, as created by \code{\link[grid]{arrow}}}
    \item{\code{force}}{Force of repulsion between overlapping text labels. Defaults
to 1.}
    \item{\code{force_pull}}{Force of attraction between a text label and its
corresponding data point. Defaults to 1.}
    \item{\code{max.time}}{Maximum number of seconds to try to resolve overlaps.
Defaults to 0.5.}
    \item{\code{max.iter}}{Maximum number of iterations to try to resolve overlaps.
Defaults to 10000.}
    \item{\code{max.overlaps}}{Exclude text labels when they overlap too many other
things. For each text label, we count how many other text labels or other
data points it overlaps, and exclude the text label if it has too many overlaps.
Defaults to 10.}
    \item{\code{nudge_x,nudge_y}}{Horizontal and vertical adjustments to nudge the
starting position of each text label. The units for \code{nudge_x} and
\code{nudge_y} are the same as for the data units on the x-axis and y-axis.}
    \item{\code{xlim,ylim}}{Limits for the x and y axes. Text labels will be constrained
to these limits. By default, text labels are constrained to the entire plot
area.}
    \item{\code{na.rm}}{If \code{FALSE} (the default), removes missing values with
a warning.  If \code{TRUE} silently removes missing values.}
    \item{\code{direction}}{"both", "x", or "y" -- direction in which to adjust position of labels}
    \item{\code{seed}}{Random seed passed to \code{\link[base]{set.seed}}. Defaults to
\code{NA}, which means that \code{set.seed} will not be called.}
    \item{\code{verbose}}{If \code{TRUE}, some diagnostics of the repel algorithm are printed}
    \item{\code{inherit.aes}}{If \code{FALSE}, overrides the default aesthetics,
rather than combining with them. This is most useful for helper functions
that define both data and aesthetics and shouldn't inherit behaviour from
the default plot specification, e.g. \code{\link[ggplot2]{borders}}.}
  }}

\item{groups}{vector of cluster labels, names contain cell names (default=NULL)}

\item{colors}{vector of numbers, which must be shown with point colors, names contain cell names (default=NULL). This argument is ignored if groups are provided.}

\item{subgroups}{subset of 'groups', selecting the cells for plot (default=NULL). Ignored if 'groups' is NULL}

\item{plot.na}{boolean/numeric Whether to plot points, for which groups / colors are missed (default=is.null(subgroups), i.e. FALSE). If plot.na passed a numeric value below 0, the NA symbols are plotted below the cells. Otherwise if values >=0, they’re plotted above the cells. Note that this argument is FALSE if 'subgroups' is NULL}

\item{min.cluster.size}{labels for all groups with number of cells fewer than this parameter are considered as missed (default=0). This argument is ignored if groups aren't provided}

\item{mark.groups}{plot cluster labels above points (default=TRUE)}

\item{show.legend}{show legend (default=FALSE)}

\item{alpha}{opacity level [0, 1] (default=0.4)}

\item{size}{point size (default=0.8)}

\item{title}{plot title (default=NULL)}

\item{plot.theme}{theme for the plot (default=NULL)}

\item{palette}{vector or list or function (default=NULL). Accepts number of colors and return list of colors (i.e. see 'colorRampPalette') (default=NULL)}

\item{color.range}{controls range, in which colors are estimated (default="symmetric"). Pass "all" to estimate range based on all values of "colors", pass "data" to estimate it only based on colors, presented in the embedding. Alternatively you can pass vector of length 2 with (min, max) values.}

\item{font.size}{font size for cluster labels (default=c(3, 7)). It can either be single number for constant font size or pair (min, max) for font size depending on cluster size}

\item{show.ticks}{show ticks and tick labels (default=FALSE)}

\item{show.labels}{show labels (default=FALSE)}

\item{legend.position}{vector with (x, y) positions of the legend (default=NULL)}

\item{legend.title}{legend title (default=NULL)}

\item{gradient.range.quantile}{Winsorization quantile for the numeric colors and gene gradient (default=1)}

\item{raster}{boolean whether layer with the points be rasterized (default=FALSE). Setting of this argument to TRUE is useful when you need to export a plot with large number of points}

\item{raster.dpi}{dpi of the rasterized plot. (default=300). Ignored if raster == FALSE.}

\item{shuffle.colors}{shuffle colors (default=FALSE)}

\item{keep.limits}{Keep axis limits from original plot (default=!is.null(subgroups)). Useful when plotting subgroups, only meaningful it plot.na=FALSE}

\item{reduction}{Reduction to use for embedding (default=NULL)}
}
\value{
ggplot2 object

ggplot2 object
}
\description{
embeddingPlot generic
Plot embedding with provided labels / colors using ggplot2

Plot embedding with provided labels / colors using ggplot2

Plot embedding from Seurat object
}
\examples{
library(sccore)
embeddingPlot(umapEmbedding, show.ticks=TRUE, show.labels=TRUE, title="UMAP embedding")

\dontrun{
# so = example Seurat object
embeddingPlot(so, groups="seurat_clusters", reduction="umap")
}

}
