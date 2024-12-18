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
#' @import parallel
#' @import RaggedExperiment
#' @export
#'
#' @examples
#'
#' data("meth_array_450k_chr14", package = "compartmap")
#' array_compartments <- getArrayABsignal(array.data.chr14, parallel = FALSE, chr = "chr14", bootstrap = FALSE, genome = "hg19", array.type = "hm450")
getArrayABsignal <- function(obj, res = 1e6, parallel = TRUE, chr = NULL,
                             targets = NULL, preprocess = TRUE, cores = 2,
                             bootstrap = TRUE, num.bootstraps = 1000,
                             genome = c("hg19", "hg38", "mm9", "mm10"),
                             other = NULL, array.type = c("hm450", "EPIC"),
                             group = FALSE, boot.parallel = TRUE, boot.cores = 2) {
  if (length(seqinfo(rowRanges(obj))) == 0) {
    message("The SummarizedExperiment you have provided has no coordinates.")
    message("Compartment extraction will fail.")
    stop("Please provide rowRanges with genomic coordinates for the object.")
  }

  # preprocess the arrays
  if (preprocess) {
    obj <- preprocessArrays(
      obj = obj,
      genome = genome, other = other,
      array.type = array.type
    )
  }

  # critical check if imputation was _not_ done
  # to send things back to beta land
  is.beta <- ifelse(min(assays(obj)$Beta, na.rm = TRUE) < 0, FALSE, TRUE)
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
      obj = obj, targets = targets, num.bootstraps = num.bootstraps,
      assay = "array", parallel = parallel, num.cores = cores
    )
  }

  if (parallel & isFALSE(group)) {
    array.compartments <- mclapply(columns, function(s) {
      obj.sub <- obj[, s]
      message("Working on ", s)
      sort(unlist(as(lapply(chr, function(c) {
        .arrayCompartments(obj.sub, obj,
          res = res,
          chr = c, targets = targets, genome = genome,
          bootstrap = bootstrap, prior.means = prior.means,
          num.bootstraps = num.bootstraps, parallel = boot.parallel,
          cores = boot.cores, group = group, bootstrap.means = bmeans
        )
      }), "GRangesList")))
    }, mc.cores = cores, mc.preschedule = F)
  }

  if (!parallel & isFALSE(group)) {
    array.compartments <- lapply(columns, function(s) {
      obj.sub <- obj[, s]
      message("Working on ", s)
      sort(unlist(as(lapply(chr, function(c) {
        .arrayCompartments(obj.sub, obj,
          res = res,
          chr = c, targets = targets, genome = genome,
          bootstrap = bootstrap, prior.means = prior.means,
          num.bootstraps = num.bootstraps, parallel = boot.parallel,
          cores = boot.cores, group = group, bootstrap.means = bmeans
        )
      }), "GRangesList")))
    })
  }

  if (parallel & isTRUE(group)) {
    array.compartments <- sort(unlist(as(mclapply(chr, function(c) {
      .arrayCompartments(obj, obj,
        res = res,
        chr = c, targets = targets, genome = genome,
        bootstrap = bootstrap, num.bootstraps = num.bootstraps, prior.means = prior.means,
        parallel = boot.parallel, cores = boot.cores, group = group, bootstrap.means = bmeans
      )
    },
    mc.cores = cores
    ), "GRangesList")))
  }

  if (!parallel & isTRUE(group)) {
    array.compartments <- sort(unlist(as(
      lapply(chr, function(c) {
        .arrayCompartments(obj, obj,
          res = res,
          chr = c, targets = targets, genome = genome,
          bootstrap = bootstrap, num.bootstraps = num.bootstraps, prior.means = prior.means,
          parallel = boot.parallel, cores = boot.cores, group = group, bootstrap.means = bmeans
        )
      }),
      "GRangesList"
    )))
  }

  # if group-level treat a little differently
  if (group) {
    return(array.compartments)
  }
  # convert to GRangesList
  array.compartments <- as(array.compartments, "CompressedGRangesList")
  # return as a RaggedExperiment
  return(RaggedExperiment(array.compartments, colData = colData(obj)))
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
#' if (require(minfiData)) {
#'   grSet <- mapToGenome(ratioConvert(preprocessNoob(RGsetEx.sub)))
#'   preprocessArrays(grSet)
#' }
#'
#' @export
preprocessArrays <- function(obj,
                             genome = c("hg19", "hg38", "mm9", "mm10"),
                             other = NULL, array.type = c("hm450", "EPIC")) {
  # make sure the input is sane
  if (!checkAssayType(obj)) stop("Input needs to be a SummarizedExperiment")

  # what genome do we have
  genome <- match.arg(genome)

  # subset the array to open sea CpGs
  obj.opensea <- filterOpenSea(obj, genome = genome, other = other)

  # convert things to M-values
  # check the names of the assays
  if (!any(getAssayNames(obj.opensea) %in% c("Beta"))) {
    stop("The assays slot should contain 'Beta' for arrays.")
  }

  # convert to M-values if beta values given
  # this should be default but allows handling if given M-values in Beta slot
  is.beta <- ifelse(min(assays(obj)$Beta, na.rm = TRUE) < 0, FALSE, TRUE)
  if (is.beta) {
    message("Converting to squeezed M-values.")
    assays(obj.opensea)$Beta <- flogit(assays(obj.opensea)$Beta)
  }

  # impute missing values if possible
  if (any(is.na(getBeta(obj.opensea)))) {
    message("Imputing missing values.")
    obj.opensea <- imputeKNN(obj.opensea, assay = "array")
  }

  return(obj.opensea)
}




# worker function
.arrayCompartments <- function(obj, original.obj, res = 1e6, chr = NULL,
                               targets = NULL,
                               genome = c("hg19", "hg38", "mm9", "mm10"),
                               prior.means = NULL, bootstrap = TRUE,
                               num.bootstraps = 1000, parallel = FALSE,
                               cores = 2, group = FALSE, bootstrap.means = NULL) {
  # this is the main analysis function for computing compartments from arrays
  # make sure the input is sane
  if (!checkAssayType(obj)) stop("Input needs to be a SummarizedExperiment")

  # what genome do we have
  genome <- match.arg(genome)

  # set the parallel back-end core number
  if (parallel) options(mc.cores = cores)

  # update
  message("Computing compartments for ", chr)
  obj <- keepSeqlevels(obj, chr, pruning.mode = "coarse")
  original.obj <- keepSeqlevels(original.obj, chr, pruning.mode = "coarse")

  # take care of the global means
  if (!is.null(prior.means)) {
    # this assumes that we've alread computed the global means
    pmeans <- as(prior.means, "GRanges")
    pmeans <- keepSeqlevels(pmeans, chr, pruning.mode = "coarse")
    # go back to a matrix
    prior.means <- as(pmeans, "matrix")
    colnames(prior.means) <- "globalMean"
  }

  # get the shrunken bins
  obj.bins <- shrinkBins(obj, original.obj,
    prior.means = prior.means, chr = chr,
    res = res, targets = targets, assay = "array",
    genome = genome, jse = TRUE
  )

  # compute correlations
  if (group) obj.cor <- getCorMatrix(obj.bins, squeeze = FALSE)
  if (isFALSE(group)) obj.cor <- getCorMatrix(obj.bins, squeeze = TRUE)
  if (any(is.na(obj.cor$binmat.cor))) {
    obj.cor$gr$pc <- matrix(rep(NA, nrow(obj.cor$binmat.cor)))
    obj.svd <- obj.cor$gr
  } else {
    # compute SVD of correlation matrix
    obj.svd <- getABSignal(obj.cor, assay = "array")
  }

  if (isFALSE(bootstrap)) {
    return(obj.svd)
  }

  # bootstrap the estimates
  # always compute confidence intervals too
  # take care of the global means
  if (bootstrap) {
    # this assumes that we've alread computed the global means
    bmeans <- as(bootstrap.means, "GRanges")
    bmeans <- keepSeqlevels(bmeans, chr, pruning.mode = "coarse")
    # go back to a matrix
    bmeans <- as(bmeans, "matrix")
    colnames(bmeans) <- rep("globalMean", ncol(bmeans))
  }

  obj.bootstrap <- bootstrapCompartments(obj,
    original.obj,
    bootstrap.samples = num.bootstraps,
    chr = chr,
    assay = "array",
    parallel = parallel,
    cores = cores,
    targets = targets,
    res = res,
    genome = genome,
    q = 0.95,
    svd = obj.svd,
    group = group,
    bootstrap.means = bmeans
  )

  # combine and return
  return(obj.bootstrap)
}
