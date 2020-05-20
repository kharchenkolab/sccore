#' @description Wrapper function around setNames 
#' @param x an object for which names attribute will be meaningful 
#' @return An object with names assigned equal to values
#' @export
sn <- function(x) {setNames(x, x)}

#' @description Extend matrix
#' @param mtx Matrix
#' @param col.names Columns that should be included in matrix
#' @return A matrix
#' @export
extendMatrix <- function(mtx, col.names) {
  new.names <- setdiff(col.names, colnames(mtx))
  ext.mtx <- matrix(0, nrow=nrow(mtx), ncol=length(new.names))
  colnames(ext.mtx) <- new.names
  return(cbind(mtx, ext.mtx)[, col.names])
}

#' @description Merge count matrices
#' @param cms List of count matrices
#' @param transposed Indicate whether cms are transposed, e.g. cells in rows and genes in columns (default=FALSE)
#' @param ... Parameters for 'plapply' function
#' @return A matrix
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