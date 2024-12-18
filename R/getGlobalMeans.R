#' Get the global means of a matrix
#'
#' @name getGlobalMeans
#'
#' @param obj Input SummarizedExperiment object
#' @param targets Column names or indices to indicate which samples to shrink towards
#' @param assay What type of assay the data are from
#'
#' @return A vector of global or targeted means
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' data("k562_scrna_chr14", package = "compartmap")
#' scrna.global.means <- getGlobalMeans(k562_scrna_chr14, assay = "rna")
#'
getGlobalMeans <- function(obj, targets = NULL,
                           assay = c("atac", "rna", "array")) {
  # match the assay arg
  assay <- match.arg(assay)

  # check the names of the assays
  if (!any(getAssayNames(obj) %in% c("counts", "Beta"))) {
    message("The assay slot should contain 'Beta' for arrays.")
    stop("The assay slot should contain 'counts' for atac/rna.")
  }

  # get the global means
  # check if shrinkage targets are being used
  if (!is.null(targets)) {
    stargets <- getShrinkageTargets(obj, targets)
    message("Using ", paste(shQuote(targets), collapse = ", "), " as shrinkage targets...")
    globalMean <- switch(assay,
      atac = matrix(rowMeans(assays(stargets)$counts, na.rm = TRUE), ncol = 1),
      rna = matrix(rowMeans(assays(stargets)$counts, na.rm = TRUE), ncol = 1),
      array = matrix(rowMeans(flogit(assays(stargets)$Beta), na.rm = TRUE), ncol = 1)
    )
  } else {
    globalMean <- switch(assay,
      atac = matrix(rowMeans(assays(obj)$counts, na.rm = TRUE), ncol = 1),
      rna = matrix(rowMeans(assays(obj)$counts, na.rm = TRUE), ncol = 1),
      array = matrix(rowMeans(flogit(assays(obj)$Beta), na.rm = TRUE), ncol = 1)
    )
  }
  colnames(globalMean) <- "globalMean"
  # coercion to get the rownames to be the GRanges coordinates
  # allows to flip back and forth between a matrix and GRanges
  # for subsetting, etc.
  rownames(globalMean) <- as.character(granges(obj))
  return(globalMean)
}

#' Pre-compute the global means for bootstrapping compartments
#'
#' @name precomputeBootstrapMeans
#'
#' @param obj Input SummarizedExperiment object
#' @param targets Optional targets to shrink towards
#' @param num.bootstraps The number of bootstraps to compute
#' @param assay What type of assay the data are from
#' @param parallel Whether to run in parallel
#' @param num.cores How many cores to use for parallel processing
#'
#' @return A matrix of bootstrapped global means
#'
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' data("k562_scrna_chr14", package = "compartmap")
#' scrna.bootstrap.global.means <- precomputeBootstrapMeans(k562_scrna_chr14, assay = "rna", num.bootstraps = 2)
#'
precomputeBootstrapMeans <- function(obj, targets = NULL, num.bootstraps = 100,
                                     assay = c("atac", "rna", "array"),
                                     parallel = FALSE, num.cores = 1) {
  # this function precomputes the bootstrapped global means
  # as a default we will make 100 bootstraps
  if (parallel) {
    bootMean <- mclapply(1:num.bootstraps, function(b) {
      message("Working on bootstrap ", b)
      if (!is.null(targets)) {
        if (length(targets) < 5) stop("Need more than 5 samples for bootstrapping to work.")
        obj <- getShrinkageTargets(obj, targets)
      }
      resamp.mat <- switch(assay,
        atac = .resampleMatrix(assays(obj)$counts),
        rna = .resampleMatrix(assays(obj)$counts),
        array = .resampleMatrix(assays(obj)$Beta)
      )
      # turn back into SummarizedExperiment
      resamp.se <- switch(assay,
        atac = SummarizedExperiment(
          assays = SimpleList(counts = resamp.mat),
          rowRanges = rowRanges(obj)
        ),
        rna = SummarizedExperiment(
          assays = SimpleList(counts = resamp.mat),
          rowRanges = rowRanges(obj)
        ),
        array = SummarizedExperiment(
          assays = SimpleList(Beta = resamp.mat),
          rowRanges = rowRanges(obj)
        )
      )
      # make sure targets is NULL since we already subset to them!
      return(getGlobalMeans(obj = resamp.se, targets = NULL, assay = assay))
    }, mc.cores = num.cores)
  } else {
    bootMean <- lapply(1:num.bootstraps, function(b) {
      message("Working on bootstrap ", b)
      if (!is.null(targets)) {
        if (length(targets) < 5) stop("Need more than 5 samples for targeted bootstrapping to be useful.")
        obj <- getShrinkageTargets(obj, targets)
      }
      resamp.mat <- switch(assay,
        atac = .resampleMatrix(assays(obj)$counts),
        rna = .resampleMatrix(assays(obj)$counts),
        array = .resampleMatrix(assays(obj)$Beta)
      )
      # turn back into SummarizedExperiment
      resamp.se <- switch(assay,
        atac = SummarizedExperiment(
          assays = SimpleList(counts = resamp.mat),
          rowRanges = rowRanges(obj)
        ),
        rna = SummarizedExperiment(
          assays = SimpleList(counts = resamp.mat),
          rowRanges = rowRanges(obj)
        ),
        array = SummarizedExperiment(
          assays = SimpleList(Beta = resamp.mat),
          rowRanges = rowRanges(obj)
        )
      )
      # make sure targets is NULL since we already subset to them!
      return(getGlobalMeans(obj = resamp.se, targets = NULL, assay = assay))
    })
  }
  return(do.call("cbind", bootMean))
}
