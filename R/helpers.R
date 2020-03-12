#' @importFrom magrittr %>% %<>% %$%
NULL

#' Parallel lapply
#' @description Parallel, optionally verbose lapply. See ?parallel::mclapply for more info.
#' @param n.cores Sumber of cores to use (default=1)
#' @param progress Show progress bar (default=F)
#' @export
plapply <- function(..., n.cores=1, progress=F, mc.preschedule=TRUE, mc.allow.recursive=TRUE) {
  if (progress && requireNamespace("pbapply", quietly=TRUE))
    return(pbapply::pblapply(..., cl=n.cores))

  if ((n.cores == 1) || !requireNamespace("parallel", quietly=TRUE))
    return(lapply(...))

  return(parallel::mclapply(..., mc.cores=n.cores, mc.preschedule=mc.preschedule, mc.allow.recursive=mc.allow.recursive))
}
