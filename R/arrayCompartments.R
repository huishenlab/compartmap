#' @title Estimate A/B compartments from methylation array data
#'
#' @description
#' \code{arrayCompartments} returns estimated A/B compartments from methylation array data.
#'
#' @param obj Input SummarizedExperiment object
#' @param res Compartment resolution in bp
#' @param chr What chromosome to work on (leave as NULL to run on all chromosomes)
#' @param group Whether to treat this as a group set of samples
#' @param targets Samples/cells to shrink towards
#' @param bootstrap Whether we should perform bootstrapping of inferred compartments
#' @param num.bootstraps How many bootstraps to run
#' @param preprocess Whether to preprocess the arrays prior to compartment inference
#' @param array.type What type of array is this ("hm450", "EPIC")
#' @param genome What genome to work on ("hg19", "hg38", "mm9", "mm10")
#' @param other Another arbitrary genome to compute compartments on
#' @param boot.parallel Whether to run the bootstrapping in parallel. See details.
#' @param BPPARAM BiocParallelParam object to use for parallelization. See details.
#'
#' @inherit scCompartments details
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
  group = FALSE,
  targets = NULL,
  bootstrap = TRUE,
  num.bootstraps = 1000,
  preprocess = TRUE,
  array.type = c("hm450", "EPIC"),
  genome = c("hg19", "hg38", "mm9", "mm10"),
  other = NULL,
  boot.parallel = TRUE,
  BPPARAM = bpparam()
) {
  verifySE(obj)
  verifyCoords(obj)
  verifyAssayNames(obj, assay = "array")
  bpparams <- get_nested_params(BPPARAM)
  check_worker_count(bpparams, group, length(chr), bootstrap)

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
    chr = chr,
    targets = targets,
    bootstrap = bootstrap,
    num.bootstraps = num.bootstraps,
    genome = genome,
    group = group,
    boot.parallel = boot.parallel,
    bpparams = bpparams
  )
}
