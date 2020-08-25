#' @importFrom magrittr %>% %<>% %$%
NULL

#' Parallel lapply
#'
#' @description Parallel, optionally verbose lapply. See ?parallel::mclapply for more info.
#' @param n.cores Number of cores to use (default=1)
#' @param progress Show progress bar (default=FALSE)
#' @examples
#' square = function(x){ x**2 }
#' plapply(1:10, square, n.cores=1, progress=TRUE)
#'
#' @return list, as returned by lapply
#' @export
plapply <- function(..., n.cores=1, progress=FALSE, mc.preschedule=TRUE, mc.allow.recursive=TRUE) {
  if (progress && requireNamespace("pbapply", quietly=TRUE)) {
    res <- pbapply::pblapply(..., cl=n.cores)
  } else if ((n.cores == 1) || !requireNamespace("parallel", quietly=TRUE)) {
  	res <- lapply(...)
  } else {
    res <- parallel::mclapply(..., mc.cores=n.cores, mc.preschedule=mc.preschedule, mc.allow.recursive=mc.allow.recursive)
  }

  is.error <- sapply(res, function(x) "try-error" %in% class(x))
  if (any(is.error)) {
    stop(paste("Errors in papply:", result[is.error]))
  }

  return(res)
}


#' Set range for values in object
#'
#' @description Changes values outside of range to min or max. Adapted from Seurat::MinMax
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


# translate multilevel segmentation into a dendrogram, with the lowest level of the dendrogram listing the cells
multi2dend <- function(cl, counts, deep=FALSE, dist='cor') {
  if(deep) {
    clf <- as.integer(cl$memberships[1,]); # take the lowest level
  } else {
    clf <- as.integer(membership(cl));
  }
  names(clf) <- names(membership(cl))
  clf.size <- unlist(tapply(clf,factor(clf,levels=seq(1,max(clf))),length))
  rowFac <- rep(NA,nrow(counts));
  rowFac[match(names(clf),rownames(counts))] <- clf;
  lvec <- colSumByFac(counts,rowFac)[-1,,drop=F];
  if(dist=='JS') {
    lvec.dist <- jsDist(t(lvec/pmax(1,Matrix::rowSums(lvec))));
  } else { # use correlation distance in log10 space
    lvec.dist <- 1-cor(t(log10(lvec/pmax(1,Matrix::rowSums(lvec))+1)))
  }
  d <- as.dendrogram(hclust(as.dist(lvec.dist),method='ward.D'))
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
  d <- dendrapply(d,addinfo,env=environment())
  attr(d,'root') <- TRUE
  return(d)
}

#' Increase resolution for a specific set of clusters
#'
#' @param con conos object
#' @param target.clusters clusters for which the resolution should be increased
#' @param clustering name of clustering in the conos object to use. Either 'clustering' or 'groups' must be provided. Default: NULL
#' @param groups set of clusters to use. Ignored if 'clustering' is not NULL. Default: NULL
#' @param method function, used to find communities. Default: leiden.community
#' @param ... additional params passed to the community function
#' @export
findSubcommunities <- function(con, target.clusters, clustering=NULL, groups=NULL, method=leiden.community, ...) {
  groups <- parseCellGroups(con, clustering, groups)

  groups.raw <- as.character(groups) %>% setNames(names(groups))
  groups <- groups[intersect(names(groups), V(con$graph)$name)]

  if(length(groups) == 0) {
    stop("'groups' not defined for graph object.")
  }

  groups <- droplevels(as.factor(groups)[groups %in% target.clusters])
  if(length(groups) == 0) {
    stop("None of 'target.clusters' can be found in 'groups'.")
  }

  subgroups <- split(names(groups), groups)
  for (n in names(subgroups)) {
    if (length(subgroups[[n]]) < 2)
      next

    new.clusts <- method(induced_subgraph(con$graph, subgroups[[n]]), ...)
    groups.raw[new.clusts$names] <- paste0(n, "_", new.clusts$membership)
  }

  return(groups.raw)
}


##' Merge into a common matrix, entering 0s for the missing ones
mergeCountMatrices <- function(cms, transposed=FALSE) {
  extendMatrix <- function(mtx, col.names) {
    new.names <- setdiff(col.names, colnames(mtx))
    ext.mtx <- Matrix::Matrix(0, nrow=nrow(mtx), ncol=length(new.names), sparse=T) %>%
      as(class(mtx)) %>% `colnames<-`(new.names)
    return(cbind(mtx, ext.mtx)[,col.names])
  }

  if (!transposed) {
    cms %<>% lapply(Matrix::t)
  }

  gene.union <- lapply(cms, colnames) %>% Reduce(union, .)
  res <- lapply(cms, extendMatrix, gene.union) %>% Reduce(rbind, .)

  if (!transposed) {
    res %<>% Matrix::t()
  }

  return(res)
}

