#' @importFrom magrittr %>% %<>% %$%
NULL

#' Parallel Lapply
#' @description parallel, optionally verbose lapply
#' @param n.cores number of cores to use
#' @param progress show progress bar
plapply <- function(..., n.cores=1, progress=F) {
  if (progress && requireNamespace("pbapply", quietly=TRUE))
    return(pbapply::pblapply(..., cl=n.cores))

  if ((n.cores == 1) || !requireNamespace("parallel", quietly=TRUE))
    return(lapply(...))

  return(parallel::mclapply(..., mc.cores=n.cores))
}
