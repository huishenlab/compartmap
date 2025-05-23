#' Invert, or "fix", compartments that have a minimum confidence score (1-min.conf)
#'
#' @name fixCompartments
#'
#' @param obj Input RaggedExperiment or output of condenseSE
#' @param min.conf Minimum confidence score to use
#' @param parallel Whether to run in parallel
#' @param cores How many cores to use if running in parallel
#'
#' @return A "fixed" set of compartments
#' @import RaggedExperiment
#' @import SummarizedExperiment
#' @export
#'
#' @examples
fixCompartments <- function(obj, min.conf = 0.8, parallel = FALSE, cores = 1) {
  # this function will invert or "fix" compartments based on bootstrapping results
  if (is(obj, "RaggedExperiment")) obj <- condenseSE(obj, sample.name = colnames(assay(obj)))

  # if we somehow only have 1 sample
  if (is(obj, "GRanges")) {
    return(flipper(obj, min.conf))
  }

  message("Fixing compartments using a minimum confidence score of ", min.conf * 100, "%")
  # go through and invert compartments based on the min.conf
  flip_compartments_lst <- mclapply(obj, flipper, min.conf, mc.cores = ifelse(parallel, cores, 1))
  names(flip_compartments_lst) <- names(obj)
  return(flip_compartments_lst)
}

#' Helper to invert, or "fix", compartments that have a minimum confidence score (1-min.conf)
#'
#' @param input_obj Input RaggedExperiment or output of condenseSE
#' @param min.conf Minimum confidence score to use
#'
#' @return A "fixed" set of compartments
#' @keywords internal
flipper <- function(input_obj, min.conf) {
  if (!any((names(mcols(input_obj)) %in% "conf.est"))) {
    stop("Bootstrapping was not performed. Cannot fix compartments.")
  }

  message("Assuming we only have a single sample to process.")
  message("Fixing compartments using a minimum confidence score of ", min.conf * 100, "%")
  invert_compartments <- apply(mcols(input_obj), 1, .inverter, min.conf)
  mcols(input_obj)$flip.compartment <- invert_compartments

  # add a new column for flipped scores
  mcols(input_obj)$flip.score <- mcols(input_obj)$score
  # flip the score
  mcols(input_obj)$flip.score[invert_compartments] <- -(mcols(input_obj)$score[invert_compartments])

  # add a new column for flipped CIs
  mcols(input_obj)$flip.conf.est <- mcols(input_obj)$conf.est
  mcols(input_obj)$flip.conf.est.upperCI <- mcols(input_obj)$conf.est.upperCI
  mcols(input_obj)$flip.conf.est.lowerCI <- mcols(input_obj)$conf.est.lowerCI

  # flip the conf.est
  mcols(input_obj)$flip.conf.est[invert_compartments] <- 1 - (mcols(input_obj)$conf.est[invert_compartments])

  # flip the upper/lowerCI
  conf.est.upperCI <- mcols(input_obj)$conf.est.upperCI[invert_compartments]
  mcols(input_obj)$flip.conf.est.upperCI[invert_compartments] <- 1 - (conf.est.upperCI)

  conf.est.lowerCI <- mcols(input_obj)$conf.est.lowerCI[invert_compartments]
  mcols(input_obj)$flip.conf.est.lowerCI[invert_compartments] <- 1 - (conf.est.lowerCI)

  return(input_obj)
}

.inverter <- function(row, min.conf) {
  row["conf.est"] < 1 - min.conf
}
