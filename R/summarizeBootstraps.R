#' Summarize the bootstrap compartment estimates and compute Agresti-Coull confidence intervals
#'
#' @name summarizeBootstraps
#'
#' @param boot.list List of bootstraps as GRanges objects
#' @param est.ab The original compartment calls
#' @param q Confidence interval to compute (0.95 for 95 percent CI)
#' @param assay Type of assay we are working with
#'
#' @return A GRanges object with bootstraps summarized
#'
#' @import SummarizedExperiment
#' @importFrom S4Vectors subjectHits queryHits
#'
#' @export
#'
#' @examples
summarizeBootstraps <- function(boot.list, est.ab, q = 0.95, assay = c("rna", "atac")) {
  # go through the estimated A/B compartments and compute proportions from the boot.list
  est.ab$score <- est.ab$pc

  is.atac_or_rna <- assay %in% c("atac", "rna")

  flog.debug("Summarizing bootstraps")

  # filter out failed compartment estimates
  boot.list <- removeEmptyBoots(boot.list)

  # summarize
  # create a dummy matrix and then rowSum them up
  boot.summary.mat.lst <- lapply(boot.list, .getSummary, est.ab = est.ab)

  # eunumerate bootstraps
  .getBootRowSums <- function(index) {
    rowSums(do.call("cbind", lapply(boot.summary.mat.lst, function(x) x[, index])))
  }
  est.ab$boot.open <- .getBootRowSums(1)
  est.ab$boot.closed <- .getBootRowSums(2)

  flog.debug("Computing Agresti-Coull 95% confidence intervals")
  .getCI(est.ab, q)
}

.getSummary <- function(gr.boot, est.ab) {
  # add the pc to the granges object
  gr.boot$score <- gr.boot$pc

  # generate a dummy GRanges object
  est.ab.dummy <- est.ab
  est.ab.dummy$boot.open <- 0
  est.ab.dummy$boot.closed <- 0

  # determine whether compartment is open and convert the boolean to 1/0 binary result for proportions
  gr.boot.isOpen <- gr.boot$compartments == "open"
  gr.boot$open <- as.integer(gr.boot.isOpen)
  gr.boot$closed <- as.integer(!gr.boot.isOpen)

  # overlap by common intervals
  ol <- findOverlaps(gr.boot, est.ab.dummy)
  mcols(est.ab.dummy)$boot.open[subjectHits(ol)] <- mcols(est.ab.dummy)$boot.open[subjectHits(ol)] + mcols(gr.boot)$open[queryHits(ol)]
  mcols(est.ab.dummy)$boot.closed[subjectHits(ol)] <- mcols(est.ab.dummy)$boot.closed[subjectHits(ol)] + mcols(gr.boot)$closed[queryHits(ol)]

  # return the dummy mcols for bootstrapped open and closed calls
  return(as.matrix(cbind(
    mcols(est.ab.dummy)$boot.open,
    mcols(est.ab.dummy)$boot.closed
  )))
}

# Add agrestiCoullCI to GRanges of compartment calls
# Uses the boot.open and boot.closed tallies across the bootstraps to get the
# successes (bootstrapped open/closed counts that match the input est.ab) and
# failures (bootstrapped open/closed counts that don't match the input est.ab)
# to pass to agrestiCoullCI
.getCI <- function(est.ab, q) {
  is.open <- est.ab$compartments == "open"
  success <- rep(NA, length(est.ab))
  failure <- rep(NA, length(est.ab))

  success[is.open] <- est.ab$boot.open[is.open]
  failure[is.open] <- est.ab$boot.closed[is.open]
  success[!is.open] <- est.ab$boot.closed[!is.open]
  failure[!is.open] <- est.ab$boot.open[!is.open]

  conf.ests <- agrestiCoullCI(success, failure, q)
  mcols(est.ab) <- cbind(mcols(est.ab), conf.ests)
  est.ab
}
