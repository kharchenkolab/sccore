#' @export
appendSpecificityMetricsToDE <- function(de.df, clusters, cluster.id, p2.counts, low.expression.threshold=0, append.auc=FALSE) {
  if (length(de.df) == 0 || nrow(de.df) == 0)
    return(de.df)

  cluster.mask <- setNames(clusters == cluster.id, names(clusters))
  if (!any(cluster.mask))
    stop("Cluster ", cluster.id, " not presented in the data")

  counts.bin <- p2.counts[names(cluster.mask), de.df$Gene, drop=F]
  counts.bin@x <- as.numeric(counts.bin@x > low.expression.threshold)
  counts.bin.sums <- Matrix::colSums(counts.bin)
  counts.bin.clust.sums <- Matrix::colSums(counts.bin * cluster.mask)

  if (append.auc) {
    if (requireNamespace("pROC", quietly = TRUE)) {
      de.df$AUC <- apply(counts.bin, 2, function(col) pROC::auc(as.integer(cluster.mask), as.integer(col),quiet=TRUE))
    } else {
      warning("You have to install pROC package to use append.auc")
    }
  }

  de.df$Specificity <- (length(cluster.mask) - counts.bin.sums) / (length(cluster.mask) - counts.bin.clust.sums)
  de.df$Precision <- counts.bin.clust.sums / counts.bin.sums
  de.df$ExpressionFraction <- Matrix::colMeans(counts.bin[cluster.mask,, drop=F])

  return(de.df)
}
