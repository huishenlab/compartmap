#' Non-parametric bootstrapping of compartments and summarization of
#' bootstraps/compute confidence intervals
#'
#' @name bootstrapCompartments
#'
#' @param obj List object of computed compartments for a sample with 'pc' and 'gr' as elements
#' @param original.obj The original, full input SummarizedExperiment of all samples/cells
#' @param BPPARAM BiocParallelParam for parallelizing bootstrapping
#' @param bootstrap.samples How many bootstraps to run
#' @param chr Which chromosome to operate on
#' @param assay What sort of assay are we working on
#' @param targets Targets to shrink towards
#' @param res The compartment resolution
#' @param genome What genome are we working on
#' @param q What sort of confidence intervals are we computing (e.g. 0.95 for 95 percentCI)
#' @param svd The original compartment calls as a GRanges object
#' @param group Whether this is for group-level inference
#' @param bootstrap.means Pre-computed bootstrap means matrix
#'
#' @return Compartment estimates with summarized bootstraps and confidence intervals
#' @import SummarizedExperiment
#'
#' @examples
#'
#' # this needs a good example
#'
#' @export
bootstrapCompartments <- function(
  obj,
  original.obj,
  BPPARAM,
  bootstrap.samples = 1000,
  chr = "chr14",
  group = FALSE,
  assay = c("rna", "atac", "array"),
  targets = NULL,
  res = 1e6,
  genome = c("hg19", "hg38", "mm9", "mm10"),
  q = 0.95,
  svd = NULL,
  bootstrap.means = NULL
) {
  # function for nonparametric bootstrap of compartments and compute 95% CIs
  # check input
  # match the assay args
  assay <- match.arg(assay)

  # if we are using targeted means
  if (!is.null(targets)) original.obj <- original.obj[, targets]

  # get the global means we are going to use
  # this could theoretically break if you ask for more bootstraps here than were pre-computed...
  # let's check for one more optimization
  if (bootstrap.samples == ncol(bootstrap.means)) {
    bmeans <- bootstrap.means
  } else {
    bmeans <- sample.int(bootstrap.means, size = bootstrap.samples, replace = FALSE)
    colnames(bmeans) <- rep("globalMean", ncol(bmeans))
  }

  # if (ncol(original.obj) < 6) stop("We need more than 5 samples to bootstrap with for the results to be meaningful.")

  # bootstrap and recompute compartments
  BiocParallel::bpprogressbar(BPPARAM) <- FALSE
  resamp.compartments <- bplapply(
    1:ncol(bmeans),
    function(b) {
      # get the shrunken bins with new global mean
      boot.mean <- as.matrix(bmeans[, b])
      colnames(boot.mean) <- "globalMean"
      s.bins <- shrinkBins(
        obj,
        original.obj,
        prior.means = boot.mean,
        chr = chr,
        res = res,
        assay = assay,
        genome = genome
      )
      cor.bins <- getCorMatrix(s.bins, squeeze = !group)

      # Stupid check for perfect correlation with global mean
      if (any(is.na(cor.bins$binmat.cor))) {
        absig <- matrix(rep(NA, nrow(cor.bins$binmat.cor)))
      } else {
        absig <- getABSignal(cor.bins, assay = assay)
      }
      return(absig)
    },
    BPPARAM = BPPARAM
  )

  # summarize the bootstraps and compute confidence intervals
  resamp.compartments <- summarizeBootstraps(resamp.compartments, svd, q = q, assay = assay)
  return(resamp.compartments)
}

# helper function to re-sample
# this was inspired by https://github.com/sgibb/bootstrap/blob/master/R/helper-functions.R
.resampleMatrix <- function(x, size = ncol(x)) {
  samp.to.select <- sample.int(ncol(x), size = size, replace = TRUE)
  return(x[, samp.to.select])
}
