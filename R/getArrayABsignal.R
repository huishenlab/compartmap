#' @title Estimate A/B compartments from methylation array data
#'
#' @description
#' \code{getArrayABsignal} returns estimated A/B compartments from methylation array data.
#'
#' @param obj Input SummarizedExperiment object
#' @param res Compartment resolution in bp
#' @param parallel Whether to run samples in parallel
#' @param chr What chromosome to work on (leave as NULL to run on all chromosomes)
#' @param targets Samples/cells to shrink towards
#' @param preprocess Whether to preprocess the arrays prior to compartment inference
#' @param cores How many cores to use when running samples in parallel
#' @param bootstrap Whether we should perform bootstrapping of inferred compartments
#' @param num.bootstraps How many bootstraps to run
#' @param genome What genome to work on ("hg19", "hg38", "mm9", "mm10")
#' @param other Another arbitrary genome to compute compartments on
#' @param array.type What type of array is this ("hm450", "EPIC")
#' @param group Whether to treat this as a group set of samples
#' @param boot.parallel Whether to run the bootstrapping in parallel
#' @param boot.cores How many cores to use for the bootstrapping
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
#'   array_compartments <- getArrayABsignal(
#'     array.data.chr14,
#'     parallel=FALSE,
#'     chr="chr14",
#'     bootstrap=FALSE,
#'     genome="hg19",
#'     array.type="hm450"
#'   )
#' }
getArrayABsignal <- function(
  obj,
  res = 1e6,
  parallel = TRUE,
  chr = NULL,
  targets = NULL,
  preprocess = TRUE,
  cores = 2,
  bootstrap = TRUE,
  num.bootstraps = 1000,
  genome = c("hg19", "hg38", "mm9", "mm10"),
  other = NULL,
  array.type = c("hm450", "EPIC"),
  group = FALSE,
  boot.parallel = TRUE,
  boot.cores = 2
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

  # gather the chromosomes we are working on
  if (is.null(chr)) {
    message("Assuming we want to process all chromosomes.")
    # get what chromosomes we want
    chr <- getChrs(obj)
  }

  # get the column names
  if (is.null(colnames(obj))) stop("colnames needs to be sample names.")
  columns <- colnames(obj)
  names(columns) <- columns

  # precompute global means
  prior.means <- getGlobalMeans(obj = obj, targets = targets, assay = "array")

  if (bootstrap) {
    message("Pre-computing the bootstrap global means.")
    bmeans <- precomputeBootstrapMeans(
      obj = obj,
      targets = targets,
      num.bootstraps = num.bootstraps,
      assay = "array",
      parallel = parallel,
      num.cores = cores
    )
  }

  if (group) {
    array.compartments.list <- mclapply(chr, function(c) {
      getCompartments(
        obj,
        obj,
        assay = "array",
        res = res,
        chr = c,
        targets = targets,
        genome = genome,
        bootstrap = bootstrap,
        num.bootstraps = num.bootstraps,
        prior.means = prior.means,
        parallel = boot.parallel,
        cores = boot.cores,
        group = group,
        bootstrap.means = bmeans
      )
    }, mc.cores = cores)
    array.compartments <- sort(unlist(as(array.compartments.list, "GRangesList")))
    return(array.compartments)
  }

  array.compartments <- mclapply(columns, function(s) {
    obj.sub <- obj[, s]
    message("Working on ", s)
    array.compartments.list <- lapply(chr, function(c) {
      getCompartments(
        obj.sub,
        obj,
        assay = "array",
        res = res,
        chr = c,
        targets = targets,
        genome = genome,
        bootstrap = bootstrap,
        prior.means = prior.means,
        num.bootstraps = num.bootstraps,
        parallel = boot.parallel,
        cores = boot.cores,
        group = group,
        bootstrap.means = bmeans
      )
    })
    sort(unlist(as(array.compartments.list, "GRangesList")))
  }, mc.cores = ifelse(parallel, cores, 1), mc.preschedule = F)

  array.compartments <- as(array.compartments, "CompressedGRangesList")
  RaggedExperiment(array.compartments, colData = colData(obj))
}

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
    message("Converting to squeezed M-values.")
    assays(obj.opensea)$Beta <- flogit(assays(obj.opensea)$Beta)
  }

  # impute missing values if possible
  if (any(is.na(minfi::getBeta(obj.opensea)))) {
    message("Imputing missing values.")
    obj.opensea <- imputeKNN(obj.opensea, assay = "array")
  }

  obj.opensea
}
