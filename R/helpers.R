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
#' @export
plapply <- function(..., n.cores=1, progress=FALSE, mc.preschedule=TRUE, mc.allow.recursive=TRUE) {
  if (progress && requireNamespace("pbapply", quietly=TRUE)) {
    return(pbapply::pblapply(..., cl=n.cores))
  }

  if ((n.cores == 1) || !requireNamespace("parallel", quietly=TRUE)) {
  	return(lapply(...))
  }

  return(parallel::mclapply(..., mc.cores=n.cores, mc.preschedule=mc.preschedule, mc.allow.recursive=mc.allow.recursive))
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