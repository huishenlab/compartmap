#' @title Estimate A/B compartments from methylation array data
#'
#' @description
#' \code{arrayCompartments} returns estimated A/B compartments from methylation array data.
#'
#' @param obj Input SummarizedExperiment object
#' @param res Compartment resolution in bp
#' @param chr What chromosome to work on (leave as NULL to run on all chromosomes)
#' @param targets Samples/cells to shrink towards
#' @param preprocess Whether to preprocess the arrays prior to compartment inference
#' @param parallel Whether to run samples in parallel
#' @param cores How many cores to use when running samples in parallel
#' @param bootstrap Whether we should perform bootstrapping of inferred compartments
#' @param num.bootstraps How many bootstraps to run
#' @param boot.parallel Whether to run the bootstrapping in parallel
#' @param boot.cores How many cores to use for the bootstrapping
#' @param genome What genome to work on ("hg19", "hg38", "mm9", "mm10")
#' @param other Another arbitrary genome to compute compartments on
#' @param group Whether to treat this as a group set of samples
#' @param array.type What type of array is this ("hm450", "EPIC")
#'
#' @return A RaggedExperiment of inferred compartments
#' @import SummarizedExperiment
#' @import RaggedExperiment
#' @importFrom parallel mclapply
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom methods as
#' @export
#'
#' @examples
#'
#' if (requireNamespace("minfi", quietly = TRUE)) {
#'   data("array_data_chr14", package = "compartmap")
#'   array_compartments <- arrayCompartments(
#'     array.data.chr14,
#'     chr="chr14",
#'     parallel=FALSE,
#'     bootstrap=FALSE,
#'     genome="hg19",
#'     array.type="hm450"
#'   )
#' }
arrayCompartments <- function(
  obj,
  res = 1e6,
  chr = NULL,
  targets = NULL,
  preprocess = TRUE,
  parallel = TRUE,
  cores = 2,
  bootstrap = TRUE,
  num.bootstraps = 1000,
  boot.parallel = TRUE,
  boot.cores = 2,
  genome = c("hg19", "hg38", "mm9", "mm10"),
  group = FALSE,
  other = NULL,
  array.type = c("hm450", "EPIC")
) {
  verifySE(obj)
  verifyCoords(obj)
  verifyAssayNames(obj, assay = "array")

  # preprocess the arrays
  if (preprocess) {
    obj <- preprocessArrays(
      obj = obj,
      genome = genome,
      other = other,
      array.type = array.type
    )
  }

  # critical check if imputation was _not_ done
  # to send things back to beta land
  is.beta <- min(assays(obj)$Beta, na.rm = TRUE) > 0
  if (!is.beta) {
    # send things back to beta land
    assays(obj)$Beta <- fexpit(assays(obj)$Beta)
  }

  getCompartments(
    obj = obj,
    assay = "array",
    res = res,
    parallel = parallel,
    chr = chr,
    targets = targets,
    cores = cores,
    bootstrap = bootstrap,
    num.bootstraps = num.bootstraps,
    boot.parallel = boot.parallel,
    boot.cores = boot.cores,
    genome = genome,
    group = group
  )
}
