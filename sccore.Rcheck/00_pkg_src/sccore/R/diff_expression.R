#' Append specificity metrics to DE
#'
#' @param de.df data.frame of differential expression values
#' @param clusters factor of clusters
#' @param cluster.id names of 'clusters' factor. If a cluster.id doesn't exist in cluster names, an error is thrown.
#' @param p2.counts counts from Pagoda2, refer to <https://github.com/kharchenkolab/pagoda2>
#' @param low.expression.threshold numeric Threshold to remove expression values (default=0). Values under this threshold are discarded. 
#' @param append.auc boolean If TRUE, append AUC values (default=FALSE)
#' @return data.frame of differential expression values with metrics attached
#' @export
appendSpecificityMetricsToDE <- function(de.df, clusters, cluster.id, p2.counts, low.expression.threshold=0, append.auc=FALSE) {
  
  if (length(de.df) == 0 || nrow(de.df) == 0) {
    return(de.df)
  }

  cluster.mask <- stats::setNames(clusters == cluster.id, names(clusters))

  if (!any(cluster.mask)) {
    stop("Cluster ", cluster.id, " not presented in the data")
  }

  counts <- p2.counts[names(cluster.mask), de.df$Gene, drop=FALSE]
  counts.bin <- counts
  counts.bin@x <- as.numeric(counts.bin@x > low.expression.threshold)
  counts.bin.sums <- Matrix::colSums(counts.bin)
  counts.bin.clust.sums <- Matrix::colSums(counts.bin * cluster.mask)

  if (append.auc) {
    de.df$AUC <- apply(counts.bin, 2, function(col) pROC::auc(as.integer(cluster.mask), as.integer(col), quiet=TRUE))
  }

  de.df$Specificity <- (length(cluster.mask) - counts.bin.sums) / (length(cluster.mask) - counts.bin.clust.sums)
  de.df$Precision <- counts.bin.clust.sums / counts.bin.sums
  de.df$ExpressionFraction <- Matrix::colMeans(counts.bin[cluster.mask,, drop=FALSE])

  return(de.df)
}
