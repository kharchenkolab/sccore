#' @importFrom magrittr %>% %<>% %$%
NULL

#' Parallel Lapply
#' @description parallel, optionally verbose lapply
#' @param n.cores number of cores to use
#' @param progress show progress bar
plapply <- function(..., n.cores=1, progress=F, mc.preschedule=TRUE, mc.allow.recursive=TRUE) {
  if (progress && requireNamespace("pbapply", quietly=TRUE))
    return(pbapply::pblapply(..., cl=n.cores))

  if ((n.cores == 1) || !requireNamespace("parallel", quietly=TRUE))
    return(lapply(...))

  return(parallel::mclapply(..., mc.cores=n.cores, mc.preschedule=mc.preschedule, mc.allow.recursive=mc.allow.recursive))
}
