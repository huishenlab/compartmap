#' Preprocess arrays for compartment inference
#'
#' @name preprocessArrays
#'
#' @param obj Input SummarizedExperiment
#' @param genome What genome are we working on ("hg19", "hg38", "mm9", "mm10")
#' @param other Another arbitrary genome to compute compartments on
#' @param array.type What type of array is this ("hm450", "EPIC")
#'
#' @return A preprocessed SummarizedExperiment to compute compartments
#' @import SummarizedExperiment
#'
#' @examples
#' if (requireNamespace("minfiData", quietly = TRUE)) {
#'   grSet <- minfi::preprocessNoob(minfiData::RGsetEx.sub) |>
#'     minfi::ratioConvert() |>
#'     minfi::mapToGenome()
#'   preprocessArrays(grSet)
#' }
#'
#' @export
preprocessArrays <- function(
  obj,
  genome = c("hg19", "hg38", "mm9", "mm10"),
  other = NULL,
  array.type = c("hm450", "EPIC")
) {
  if (!requireNamespace("minfi", quietly = TRUE)) {
    stop("The minfi package must be installed for this functionality")
  }

  # what genome do we have
  genome <- match.arg(genome)

  # subset the array to open sea CpGs
  obj.opensea <- filterOpenSea(obj, genome = genome, other = other)
  verifyAssayNames(obj.opensea, assay = "array")

  # convert to M-values if beta values given
  # this should be default but allows handling if given M-values in Beta slot
  is.beta <- min(assays(obj)$Beta, na.rm = TRUE) > 0
  if (is.beta) {
    flog.debug("Converting to squeezed M-values.")
    assays(obj.opensea)$Beta <- flogit(assays(obj.opensea)$Beta)
  }

  # impute missing values if possible
  if (any(is.na(minfi::getBeta(obj.opensea)))) {
    flog.debug("Imputing missing values.")
    obj.opensea <- imputeKNN(obj.opensea, assay = "array")
  }

  obj.opensea
}
