#' @useDynLib sccore
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' @importFrom magrittr %>% %<>% %$%
#' @import igraph
#' @importFrom stats as.dendrogram as.dist is.leaf median na.omit quantile setNames
NULL

#' Parallel, optionally verbose lapply. See ?parallel::mclapply for more info.
#'
#' @param ... Additional arguments passed to mclapply(), lapply(), or pbapply::pblapply()
#' @param progress Show progress bar via pbapply (default=FALSE)
#' @param n.cores Number of cores to use (default=parallel::detectCores())
#' @param mc.preschedule See ?parallel::mclapply (default=FALSE) If TRUE then the computation is first divided to (at most) as many jobs are there are cores and then the jobs are started, each job possibly covering more than one value. If FALSE, then one job is forked for each value of X. The former is better for short computations or large number of values in X, the latter is better for jobs that have high variance of completion time and not too many values of X compared to mc.cores.
#' @return list, as returned by lapply
#' @examples
#' square = function(x){ x**2 }
#' plapply(1:10, square, n.cores=1, progress=TRUE)
#'
#' @export
plapply <- function(..., progress=FALSE, n.cores=parallel::detectCores(), mc.preschedule=FALSE) {
  if (progress && requireNamespace("pbapply", quietly=TRUE)) {
    result <- pbapply::pblapply(..., cl=n.cores)
  } else if(n.cores>1) {
    result <- parallel::mclapply(..., mc.cores=n.cores, mc.preschedule=mc.preschedule)
  } else {
    # fall back on lapply
    result <- lapply(...)
  }

  is.error <- (sapply(result, class) == "try-error")
  if (any(is.error)) {
    stop(paste("Errors in plapply:", result[is.error]))
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
setMinMax <- function(obj, min, max) {
  obj[obj<min] <- min
  obj[obj>max] <- max
  return(obj)
}


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
  lvec <- colSumByFac(counts, rowFac)[-1,, drop=FALSE]
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


#' Increase resolution for a specific set of clusters
#'
#' @param con conos object, from <https://github.com/kharchenkolab/conos>, "Joint analysis of heterogeneous single-cell RNA-seq dataset collections", DOI: 10.1038/s41592-019-0466-z
#' @param target.clusters Clusters for which the resolution should be increased
#' @param clustering Name of clustering in the conos object to use (default=NULL). Either 'clustering' or 'groups' must be provided.
#' @param groups Set of clusters to use (default=NULL). Ignored if 'clustering' is not NULL.
#' @param method Function, used to find communities (default=igraph::cluster_louvain)
#' @param ... Additional params passed to the community function
#' @return The input conos object "subset" with the input target.clusters, thereby resulting in a magnified view of these clusters
#' @export
findSubcommunities <- function(con, target.clusters, clustering=NULL, groups=NULL, method=igraph::cluster_louvain, ...) {

  parseCellGroups <- function(con, clustering, groups) {

    if (!is.null(groups)) {
      if (!any(names(groups) %in% names(con$getDatasetPerCell()))){
        stop("'groups' aren't defined for any of the cells.")
      }
      
      return(groups)
    }

    if (is.null(clustering)) {
      if (length(con$clusters) > 0){
        return(con$clusters[[1]]$groups)
      }

      stop("Either 'groups' must be provided or the conos object must have some clustering estimated")
    }
    if(is.null(clusters[[clustering]])){
      stop(paste("clustering",clustering,"doesn't exist, run findCommunity() first"))
    }

    return(con$clusters[[clustering]]$groups)
  }

  groups <- parseCellGroups(con=con, clustering=clustering, groups=groups)

  groups.raw <- as.character(groups) %>% stats::setNames(names(groups))
  groups <- groups[intersect(names(groups), V(con$graph)$name)]

  if (length(groups) == 0) {
    stop("'groups' not defined for graph object.")
  }

  groups <- droplevels(as.factor(groups)[groups %in% target.clusters])
  if (length(groups) == 0) {
    stop("None of 'target.clusters' can be found in 'groups'.")
  }

  subgroups <- split(names(groups), groups)
  for (n in names(subgroups)) {
    if (length(subgroups[[n]]) < 2){
      next
    }

    new.clusts <- method(induced_subgraph(con$graph, subgroups[[n]]), ...)
    groups.raw[new.clusts$names] <- paste0(n, "_", new.clusts$membership)
  }

  return(groups.raw)
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
  stats::setNames(x, x)
}


#' Extend matrix to include new columns in matrix
#'
#' @param mtx Matrix
#' @param col.names Columns that should be included in matrix
#' @return Matrix with new columns but rows retained
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


