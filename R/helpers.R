#' @importFrom magrittr %>% %<>% %$%
NULL


#' Parallel lapply
#'
#' @description Parallel, optionally verbose lapply. See ?parallel::mclapply for more info.
#' @param progress Show progress bar via pbapply (default=FALSE)
#' @param n.cores Number of cores to use (default=1)
#' @param mc.preschedule See ?parllel::mclapply (default=FALSE) If TRUE then the computation is first divided to (at most) as many jobs are there are cores and then the jobs are started, each job possibly covering more than one value. If FALSE, then one job is forked for each value of X. The former is better for short computations or large number of values in X, the latter is better for jobs that have high variance of completion time and not too many values of X compared to mc.cores.
#' @examples
#' square = function(x){ x**2 }
#' plapply(1:10, square, n.cores=1, progress=TRUE)
#'
#' @return list, as returned by lapply
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
