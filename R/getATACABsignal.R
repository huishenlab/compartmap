#' @title Estimate A/B compartments from ATAC-seq data
#'
#' @description
#' \code{getATACABsignal} returns estimated A/B compartments from ATAC-seq data.
#'
#' @param obj Input SummarizedExperiment object
#' @param res Compartment resolution in bp
#' @param parallel Whether to run samples in parallel
#' @param chr What chromosome to work on (leave as NULL to run on all chromosomes)
#' @param targets Samples/cells to shrink towards
#' @param cores How many cores to use when running samples in parallel
#' @param bootstrap Whether we should perform bootstrapping of inferred compartments
#' @param num.bootstraps How many bootstraps to run
#' @param genome What genome to work on ("hg19", "hg38", "mm9", "mm10")
#' @param other Another arbitrary genome to compute compartments on
#' @param group Whether to treat this as a group set of samples
#' @param boot.parallel Whether to run the bootstrapping in parallel
#' @param boot.cores How many cores to use for the bootstrapping
#'
#' @return A RaggedExperiment of inferred compartments
#' @import SummarizedExperiment
#' @import RaggedExperiment
#' @importFrom parallel mclapply
#' @importFrom methods as
#' @export
#'
#' @aliases getRNAABsignal
#'
#' @examples
#' if (requireNamespace("csaw", quietly = TRUE)) {
#'   data("k562_scatac_chr14", package = "compartmap")
#'   atac_compartments <- getATACABsignal(
#'     k562_scatac_chr14,
#'     parallel = FALSE,
#'     chr = "chr14",
#'     bootstrap = FALSE,
#'     genome = "hg19",
#'     group = TRUE
#'   )
#' }
getATACABsignal <- function(
  obj,
  res = 1e6,
  parallel = FALSE,
  chr = NULL,
  targets = NULL,
  cores = 2,
  bootstrap = TRUE,
  num.bootstraps = 100,
  genome = c("hg19", "hg38", "mm9", "mm10"),
  other = NULL,
  group = FALSE,
  boot.parallel = FALSE,
  boot.cores = 2
) {
  chr <- chr %||%
    {
      message("Processing all chromosomes")
      getChrs(obj)
    }

  if (is.null(colnames(obj))) stop("colnames needs to be sample names.")
  columns <- colnames(obj)
  names(columns) <- columns

  prior.means <- getGlobalMeans(obj = obj, targets = targets, assay = "atac")

  if (bootstrap) {
    message("Pre-computing the bootstrap global means.")
    bmeans <- precomputeBootstrapMeans(
      obj = obj,
      targets = targets,
      num.bootstraps = num.bootstraps,
      assay = "atac",
      parallel = parallel,
      num.cores = cores
    )
  }

  if (group) {
    atac.compartments.list <- mclapply(chr, function(c) {
      getCompartments(
        obj,
        obj,
        assay = "atac",
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
    }, mc.cores = ifelse(parallel, cores, 1))

    atac.compartments <- sort(unlist(as(atac.compartments.list, "GRangesList")))
    return(atac.compartments)
  }

  atac.compartments <- mclapply(columns, function(s) {
    obj.sub <- obj[, s]

    message("Working on ", s)
    atac.compartments.list <- lapply(chr, function(c) {
      getCompartments(
        obj.sub,
        obj,
        assay = "atac",
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

    sort(unlist(as(atac.compartments.list, "GRangesList")))
  }, mc.cores = ifelse(parallel, cores, 1))

  atac.compartments <- as(atac.compartments, "CompressedGRangesList")
  RaggedExperiment(atac.compartments, colData = colData(obj))
}

#' @describeIn getATACABsignal Alias for getATACABsignal
#'
getRNAABsignal <- getATACABsignal
