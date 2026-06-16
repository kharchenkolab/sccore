## Benjamini-Hochberg adjustment, optionally in log space (so very small p-values keep their precision).
## Ported from pagoda2 so the shared marker core does not depend on pagoda2 internals.
bh.adjust <- function(x, log = FALSE) {
  nai <- which(!is.na(x))
  ox <- x
  x <- x[nai]
  id <- order(x, decreasing = FALSE)
  if (log) {
    q <- x[id] + log(length(x) / seq_along(x))
  } else {
    q <- x[id] * length(x) / seq_along(x)
  }
  a <- rev(cummin(rev(q)))[order(id)]
  ox[nai] <- a
  ox
}

#' Differentially expressed genes per group, from a count matrix
#'
#' The shared core of marker detection: for each group it runs a Mann-Whitney (Wilcoxon rank-sum) test of
#' that group against all other cells, gene by gene, on an in-memory (cells x genes) expression matrix, and
#' reports a standardized, multiple-testing-adjusted Z score together with the log2 fold change. This is the
#' computation underlying \code{pagoda2}'s \code{getDifferentialGenes()} and \code{conos}'s
#' \code{runMarkers()}; keeping it here lets both call the same code on whatever sample type they hold
#' (pagoda2, Seurat, ...). Disk-backed containers should materialize the block themselves (e.g. via a
#' streaming accessor) and pass the resulting matrix in.
#'
#' @param cm a cells x genes expression matrix (sparse \code{dgCMatrix}); rows are cells, columns are genes.
#' @param groups a factor (named by cell, or aligned to \code{rownames(cm)}) assigning each cell to a group;
#'   \code{NA} cells are dropped from the comparison.
#' @param z.threshold numeric keep only genes whose (signed or absolute) Z exceeds this (default=0).
#' @param upregulated.only boolean if TRUE keep only up-regulated genes (Z > threshold); else |Z| (default=FALSE).
#' @param append.specificity.metrics boolean append Specificity/Precision/ExpressionFraction (and AUC) via
#'   \code{appendSpecificityMetricsToDE()} (default=FALSE).
#' @param append.auc boolean if appending specificity metrics, also append AUC (default=FALSE).
#' @param n.cores integer number of cores for the specificity-metrics step (default=1).
#' @return a named list (one element per group level) of data.frames with columns \code{Z}, \code{M}
#'   (log2 fold change), \code{highest} (is this the group of maximal fold change for the gene), \code{fe}
#'   (fraction of group cells expressing) and \code{Gene}, ordered by decreasing Z.
#' @export
matrixDE <- function(cm, groups, z.threshold = 0, upregulated.only = FALSE,
                     append.specificity.metrics = FALSE, append.auc = FALSE, n.cores = 1) {
  if (!inherits(cm, "dgCMatrix")) {
    cm <- methods::as(methods::as(cm, "CsparseMatrix"), "dgCMatrix")
  }
  cols <- groups
  if (!is.null(names(cols))) {
    cols <- cols[match(rownames(cm), names(cols))]
  }
  valid <- !is.na(cols)
  if (!any(valid)) stop("No cells with non-missing groups are present in the matrix")
  if (!all(valid)) cm <- cm[valid, , drop = FALSE]
  cols <- as.factor(droplevels(as.factor(cols[valid])))

  # Mann-Whitney (Wilcoxon rank-sum) of each group vs the rest, gene by gene
  xr <- sparse_matrix_column_ranks(cm)
  grs <- colSumByFactor(xr, cols)[-1, , drop = FALSE]                  # rank sums per group
  xr@x <- numeric(length(xr@x)) + 1
  gnzz <- colSumByFactor(xr, cols)[-1, , drop = FALSE]                 # non-zero counts per group
  group.size <- as.numeric(tapply(cols, cols, length))[1:nrow(gnzz)]
  group.size[is.na(group.size)] <- 0
  gnz <- (group.size - gnzz)
  zero.ranks <- (nrow(xr) - diff(xr@p) + 1) / 2
  ustat <- t((t(gnz) * zero.ranks)) + grs - group.size * (group.size + 1) / 2
  n1n2 <- group.size * (nrow(cm) - group.size)
  usigma <- sqrt((nrow(cm) + 1 - (gnz^3 - gnz) / (nrow(cm) * (nrow(cm) - 1))) * n1n2 / 12)
  x <- t((ustat - n1n2 / 2) / usigma)                                  # standardized U -> Z
  x <- matrix(stats::qnorm(bh.adjust(stats::pnorm(as.numeric(abs(x)), lower.tail = FALSE, log.p = TRUE), log = TRUE),
                           lower.tail = FALSE, log.p = TRUE), ncol = ncol(x)) * sign(x)
  rownames(x) <- colnames(cm)
  colnames(x) <- levels(cols)[1:ncol(x)]

  # log2 fold change, fraction expressing, group of maximal fold change
  log.gene.av <- log2(Matrix::colMeans(cm))
  group.gene.av <- colSumByFactor(cm, cols)[-1, , drop = FALSE] / (group.size + 1)
  log2.fold.change <- log2(t(group.gene.av)) - log.gene.av
  f.expressing <- t(gnzz / group.size)
  max.group <- max.col(log2.fold.change)

  ds <- lapply(1:ncol(x), function(i) {
    z <- x[, i]
    vi <- which((if (upregulated.only) z else abs(z)) >= z.threshold)
    r <- data.frame(Z = z[vi], M = log2.fold.change[vi, i], highest = max.group[vi] == i,
                    fe = f.expressing[vi, i], Gene = rownames(x)[vi], stringsAsFactors = FALSE)
    rownames(r) <- r$Gene
    r[order(r$Z, decreasing = TRUE), ]
  })
  names(ds) <- colnames(x)

  if (append.specificity.metrics) {
    nm <- stats::setNames(names(ds), names(ds))
    ds <- plapply(nm, function(n) appendSpecificityMetricsToDE(ds[[n]], cols, n, p2.counts = cm, append.auc = append.auc),
                  n.cores = n.cores, progress = FALSE)
  }
  ds
}

#' Append specificity metrics to DE
#'
#' @param de.df data.frame of differential expression values
#' @param clusters factor of clusters
#' @param cluster.id names of 'clusters' factor. If a cluster.id doesn't exist in cluster names, an error is thrown.
#' @param p2.counts counts from Pagoda2, refer to <https://github.com/kharchenkolab/pagoda2>
#' @param low.expression.threshold numeric Threshold to remove expression values (default=0). Values under this threshold are discarded.
#' @param append.auc boolean If TRUE, append AUC values (default=FALSE)
#' @return data.frame of differential expression values with metrics attached
#'
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
    ## AUC of each (binary) expression indicator vs cluster membership, in closed form via the
    ## Mann-Whitney identity -- replaces a per-column pROC::auc() loop. counts.bin is binary, so with
    ## a = #expressed in-cluster, c = #expressed out-of-cluster, and b/d the non-expressed complements,
    ## AUC = (a*d + 0.5*(a*c + b*d)) / (n.pos*n.neg). This is the fixed up-regulation direction (higher
    ## expression -> in-cluster; equals pROC::auc(direction="<")): >0.5 = up-marker, <0.5 = down,
    ## consistent with the other (presence-oriented) metrics here. Constant predictors give 0.5.
    n.pos <- sum(cluster.mask)
    n.neg <- length(cluster.mask) - n.pos
    a <- counts.bin.clust.sums
    c <- counts.bin.sums - a
    b <- n.pos - a
    d <- n.neg - c
    de.df$AUC <- (a * d + 0.5 * (a * c + b * d)) / (n.pos * n.neg)
  }

  de.df$Specificity <- (length(cluster.mask) - counts.bin.sums) / (length(cluster.mask) - counts.bin.clust.sums)
  de.df$Precision <- counts.bin.clust.sums / counts.bin.sums
  de.df$ExpressionFraction <- Matrix::colMeans(counts.bin[cluster.mask,, drop=FALSE])

  return(de.df)
}

#' Collapse count matrices by cell type, given min/max number of cells
#'
#' @param cm count matrix
#' @param groups factor specifying cell types
#' @param min.cell.count numeric Minimum number of cells to include (default=10)
#' @param max.cell.count numeric Maximum number of cells to include (default=Inf). If Inf, there is no maximum.
#' @return Subsetted factor of collapsed cells by type, with NA cells omitted
#' @export
collapseCellsByType <- function(cm, groups, min.cell.count=10, max.cell.count=Inf) {
  groups <- as.factor(groups);
  cl <- setNames(factor(groups[match(rownames(cm),names(groups))],levels=levels(groups)),rownames(cm));
  if(is.finite(max.cell.count)) {
    vc <- unlist(tapply(names(cl),cl,function(nn) {
      if(length(nn) > max.cell.count) sample(nn,max.cell.count) else nn
    }))
    cl <- cl[names(cl) %in% vc]
    cm <- cm[names(cl),]
  }

  tc <- colSumByFactor(cm,cl);
  tc <- tc[-1,,drop=FALSE]  # omit NA cells
  tc[table(cl)>=min.cell.count,,drop=FALSE]
}

