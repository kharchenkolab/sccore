#' @useDynLib sccore
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' @importFrom magrittr %>% %<>% %$%
#' @import igraph
#' @importFrom stats as.dendrogram as.dist is.leaf median na.omit quantile setNames
NULL

#' Parallel, optionally verbose lapply. See ?parallel::mclapply for more info.
#'
#' @param ... Additional arguments passed to mclapply(), lapply(), or pbmcapply::pbmclapply()
#' @param progress Show progress bar via pbmcapply::pbmclapply() (default=FALSE).
#' @param fail.on.error boolean Whether to fail and report and error (using stop()) as long as any of the individual tasks has failed (default =FALSE)
#' @param n.cores Number of cores to use (default=parallel::detectCores()). When n.cores=1, regular lapply() is used. Note: doesn't work when progress=TRUE
#' @inheritParams parallel::mclapply
#' @return list, as returned by lapply
#' @examples
#' square = function(x){ x**2 }
#' plapply(1:10, square, n.cores=1, progress=TRUE)
#'
#' @export
plapply <- function(..., progress=FALSE, n.cores=parallel::detectCores(), mc.preschedule=FALSE, fail.on.error=FALSE) {
  if (progress) {
    result <- pbmcapply::pbmclapply(..., mc.cores=n.cores, mc.preschedule=mc.preschedule)
  } else if(n.cores > 1) {
    result <- parallel::mclapply(..., mc.cores=n.cores, mc.preschedule=mc.preschedule)
  } else {
    # fall back on lapply
    result <- lapply(...)
  }

  if(fail.on.error) {
    is.error <- (sapply(result, class) == "try-error")
    if (any(is.error)) {
      stop(paste("Errors in plapply:", result[is.error]))
    }
  }

  return(result)
}


#' Set range for values in object. Changes values outside of range to min or max. Adapted from Seurat::MinMax
#'
#' @param obj Object to manipulate
#' @param min Minimum value
#' @param max Maximum value
#' @return An object with the same dimensions as input but with altered range in values
#' @examples
#' example_matrix =  matrix(rep(c(1:5), 3), 5)
#' setMinMax(example_matrix, 2, 4)
#'
#' @export
setMinMax <- function(obj, min, max) { obj %>% pmax(min) %>% pmin(max) }



#' Translate multilevel segmentation into a dendrogram, with the lowest level of the dendrogram listing the cells
#'
#' @param cl igraph communities object, returned from igraph community detection functions
#' @param counts dgCmatrix of counts
#' @param deep boolean If TRUE, take (cl$memberships[1,]). Otherwise, uses as.integer(membership(cl)) (default=FALSE)
#' @param dist Distance metric used (default='cor'). Eiether 'cor' for the correlation distance in log10 space, or 'JS' for the Jensen–Shannon distance metric (i.e. the square root of the Jensen–Shannon divergence)
#' @return resulting dendrogram
#' @export
multi2dend <- function(cl, counts, deep=FALSE, dist='cor') {
  if(deep) {
    clf <- as.integer(cl$memberships[1,]); # take the lowest level
  } else {
    clf <- as.integer(membership(cl));
  }
  names(clf) <- names(membership(cl))
  clf.size <- unlist(tapply(clf, factor(clf, levels=seq(1,max(clf))), length))
  rowFac <- rep(NA,nrow(counts));
  rowFac[match(names(clf),rownames(counts))] <- clf;
  lvec <- colSumByFactor(counts, rowFac)[-1,, drop=FALSE]
  if(dist=='JS') {
    lvec.dist <- jsDist(t(lvec/pmax(1, Matrix::rowSums(lvec))));
  } else { # use correlation distance in log10 space
    lvec.dist <- 1-stats::cor(t(log10(lvec/pmax(1,Matrix::rowSums(lvec))+1)))
  }
  d <- as.dendrogram(stats::hclust(as.dist(lvec.dist),method='ward.D'))
  # add cell info to the laves
  addinfo <- function(l,env) {
    v <- as.integer(mget("index",envir=env,ifnotfound=0)[[1]])+1;
    attr(l,'nodeId') <- v
    assign("index",v,envir=env)
    attr(l,'nCells') <- sum(clf.size[as.integer(unlist(l))]);
    if(is.leaf(l)) {
      attr(l,'cells') <- names(clf)[clf==attr(l,'label')];
    }
    attr(l,'root') <- FALSE
    return(l);
  }
  d <- stats::dendrapply(d,addinfo,env=environment())
  attr(d,'root') <- TRUE
  return(d)
}

#' Set names equal to values, a stats::setNames wrapper function
#'
#' @param x an object for which names attribute will be meaningful
#' @return An object with names assigned equal to values
#' @examples
#' vec = c(1, 2, 3, 4)
#' sn(vec)
#'
#' @export
sn <- function(x) {
  names(x) <- x
  return(x)
}

#' Extend matrix to include new columns in matrix
#'
#' @param mtx Matrix
#' @param col.names Columns that should be included in matrix
#' @return Matrix with new columns but rows retained
#' @examples
#' library(dplyr)
#' geneUnion <- lapply(conosClusterList, colnames) %>% Reduce(union, .)
#' extendMatrix(conosClusterList[[1]], col.names=geneUnion)
#'
#' @export
extendMatrix <- function(mtx, col.names) {
  new.names <- setdiff(col.names, colnames(mtx))
  ## if all col.names already included in matrix, don't extend
  if (identical(new.names, character(0))){
    return(mtx)
  }
  ext.mtx <- matrix(0, nrow=nrow(mtx), ncol=length(new.names))
  colnames(ext.mtx) <- new.names
  return(cbind(mtx, ext.mtx)[, col.names])
}

#' Merge list of count matrices into a common matrix, entering 0s for the missing entries
#'
#' @param cms List of count matrices
#' @param transposed boolean Indicate whether 'cms' is transposed, e.g. cells in rows and genes in columns (default=FALSE)
#' @param ... Parameters for 'plapply' function
#' @return A merged extended matrix, with 0s for missing entries
#' @examples
#' mergeCountMatrices(conosClusterList, n.cores=1)
#' ## 12 x 67388 sparse Matrix of class "dgCMatrix"
#'
#' @export
mergeCountMatrices <- function(cms, transposed=FALSE, ...) {
  if (!transposed) {
    cms %<>% plapply(Matrix::t, ...)
  }

  gene.union <- lapply(cms, colnames) %>% Reduce(union, .)

  res <- plapply(cms, extendMatrix, gene.union, ...) %>% Reduce(rbind, .)
  if (!transposed) {
    res %<>% Matrix::t()
  }
  return(res)
}
