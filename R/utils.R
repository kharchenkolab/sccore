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
#' @return Matrix with new columns
#' @export
extendMatrix <- function(mtx, col.names) {
  new.names <- setdiff(col.names, colnames(mtx))
  ext.mtx <- matrix(0, nrow=nrow(mtx), ncol=length(new.names))
  colnames(ext.mtx) <- new.names
  return(cbind(mtx, ext.mtx)[, col.names])
}


#' Merge list of count matrices
#'
#' @param cms List of count matrices
#' @param transposed boolean Indicate whether 'cms' is transposed, e.g. cells in rows and genes in columns (default=FALSE)
#' @param ... Parameters for 'plapply' function
#' @return Mrged matrix
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
