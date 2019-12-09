#' @export
appendSpecificityMetricsToDE <- function(de.df, clusters, cluster.id, p2.counts, low.expression.threshold=0, append.auc=FALSE) {
  cluster.mask <- setNames(clusters == cluster.id, names(clusters))

  counts.bin <- (p2.counts[names(cluster.mask), de.df$Gene, drop=F] > low.expression.threshold)
  counts.bin.sums <- Matrix::colSums(counts.bin)
  counts.bin.clust.sums <- Matrix::colSums(counts.bin & cluster.mask)

  if (append.auc) {
    if (requireNamespace("pROC", quietly = TRUE)) {
      de.df$AUC <- apply(counts.bin, 2, function(col) pROC::auc(as.integer(cluster.mask), as.integer(col), quiet=T))
    } else {
      warning("You have to install pROC package to use append.auc")
    }
  }

  de.df$Specificity <- (length(cluster.mask) - counts.bin.sums) / (length(cluster.mask) - counts.bin.clust.sums)
  de.df$Precision <- counts.bin.clust.sums / counts.bin.sums
  de.df$ExpressionFraction <- Matrix::colMeans(counts.bin[cluster.mask,, drop=F])

  return(de.df)
}
