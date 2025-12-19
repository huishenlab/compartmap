#' Check and tell users whether input BPPARAM are optimal given group/sc
#' inference and bootstrapping
#' @keywords internal
check_optim <- function(workers, group, chr_count, bootstrap) {
  if (group & chr_count < workers[1]) {
    flog.info(
      "Grouped inference with more outer workers than chromosomes leaves %d of %d workers unused",
      workers[1] - chr_count,
      workers[1]
    )
    if (bootstrap) {
      flog.info("Consider using a single core for the outer worker and more cores for the inner bootstrap worker")
    }
  }
  if (!group & bootstrap & workers[1] < workers[2]) {
    flog.info("More outer (column-wise) than inner (bootstrap) workers is faster for single-cell inference")
  }
}

#' Check that the number of requested workers is valid
#' @keywords internal
check_worker_count <- function(bpparam, group, chr_count, bootstrap, avail_workers = parallelly::availableCores()) {
  workers <- get_bpnworkers(bpparam)
  total <- required_workers(workers)
  if (verify_workers(total)) {
    check_optim(workers, group, chr_count, bootstrap)
    return(TRUE)
  }

  avail_msg <- sprintf("but your system has only %d cores", avail_workers)
  info_msg <- "See parallelly::availableCores(which = 'all') for more information on available resources"
  if (workers[1] == 1 | workers[2] == 1) {
    msg <- sprintf(
      "Requested %d %s workers %s\n%s",
      max(workers),
      ifelse(workers[1] == 1, "inner", "outer"),
      avail_msg,
      info_msg
    )
  } else {
    msg <- sprintf(
      "Requested %1$d outer and %2$d inner workers that require %3$d total workers (%1$d + (%1$d x %2$d)) %4$s\n%5$s",
      workers[1],
      workers[2],
      total,
      avail_msg,
      info_msg
    )
  }
  stop(msg)
}

#' Get the total number of BiocParallelParam workers that will get used
#' From a BiocParallelParam or list of 2 BiocParallelParam objects
#' @keywords internal
get_bpnworkers <- function(bp) {
  workers <- bpnworkers.list(bp)
}

#' Return the number of workers in a list of BiocParallelParam objects
#' @param List of BiocParallelParam objects
#' @importFrom BiocParallel bpnworkers
#' @return A vector of the `bpnworkers` count in each list element
#' @keywords internal
bpnworkers.list <- function(bplist) {
  unlist(Map(bpnworkers, bplist))
}

required_workers <- function(workers) {
  if (workers[1] == 1) {
    return(workers[2])
  }
  if (workers[2] == 1) {
    return(workers[1])
  }
  sum(Reduce(`*`, workers), workers[1])
}

#' Verify that the input BiocParallelParam is valid
#' @param A BiocParallelParam or list of 2 BiocParallelParam objects
#' @importFrom BiocParallel bpnworkers
#' @return TRUE if the total `bpnworkers` in the input does not exceed
#' available resources as defined by `parallelly::availableCores()`
#' @keywords internal
verify_bp <- function(bp) {
  verify_workers(get_bpnworkers(bp))
}

#' Verify that requested thread count is not higher than available
#' @param thread_count The number of workers to check availability
#' @return TRUE if the requested `thread_count` does not exceed available
#' resources as defined by `parallelly::availableCores()`
#' @keywords internal
verify_workers <- function(n_workers) {
  avail_workers <- parallelly::availableCores()
  n_workers <= avail_workers
}

#' Set outer and inner params for nester parallelization
#' The outer param is across the input samples/columns and the second is for
#' bootstrapping. If `boot.parallel` is FALSE, the inner param is set to
#' `SerialParam`.
#' @importFrom BiocParallel SerialParam
#' @keywords internal
get_nested_params <- function(BPPARAM, boot.parallel) {
  stopifnot("Only two BiocParallelParam objects can be used" = length(BPPARAM) <= 2)
  single_param <- length(BPPARAM) == 1
  BPPARAM <- if (single_param & is.list(BPPARAM)) BPPARAM[[1]] else BPPARAM

  if (boot.parallel) {
    if (single_param) {
      return(list(outer = BPPARAM, inner = BPPARAM))
    } else {
      return(list(outer = BPPARAM[[1]], inner = BPPARAM[[2]]))
    }
  }

  if (single_param) {
    list(outer = BPPARAM, inner = SerialParam())
  } else {
    list(outer = BPPARAM[[1]], inner = SerialParam())
  }
}
