#' @title Employ an eBayes shrinkage approach for bin-level estimates for A/B inference
#'
#' @description
#' \code{shrinkBins} returns shrunken bin-level estimates
#'
#' @details This function computes shrunken bin-level estimates using a
#' James-Stein estimator (JSE), reformulated as an eBayes procedure. JSE can be
#' used only if at least 4 targets are provided - any less and `shrinkBins`
#' will fall back to using Bayes rule which will probably not be great but it
#' won't explode and may provide some reasonable results anyway
#'
#' @param x Input SummarizedExperiment object
#' @param original.x Full sample set SummarizedExperiment object
#' @param prior.means The means of the bin-level prior distribution
#' @param chr The chromosome to operate on
#' @param res Resolution to perform the binning
#' @param targets The column/sample/cell names to shrink towards
#' @param jse Whether to use a James-Stein estimator (default is TRUE)
#' @param assay What assay type this is ("rna", "atac", "array")
#' @param genome What genome are we working with ("hg19", "hg38", "mm9", "mm10")
#'
#' @return A list object to pass to getCorMatrix
#'
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @importFrom stats median sd
#'
#' @export
#'
#' @examples
#' data("k562_scrna_chr14", package = "compartmap")
#' shrunken.bin.scrna <- shrinkBins(
#'   x = k562_scrna_chr14,
#'   original.x = k562_scrna_chr14,
#'   chr = "chr14", assay = "rna"
#' )
#'
shrinkBins <- function(
  x,
  original.x,
  prior.means = NULL,
  chr = NULL,
  res = 1e6,
  targets = NULL,
  jse = TRUE,
  assay = c("rna", "atac", "array"),
  genome = c("hg19", "hg38", "mm9", "mm10")
) {
  assay <- match.arg(assay)
  genome <- match.arg(genome)
  verifySE(x)

  if (!is.null(targets)) {
    target.count <- length(targets)
    if (target.count == 1) {
      stop("Cannot perform targeted bin-level shrinkage with one target sample.")
    } else if (target.count < 4) {
      message("Number of means fewer than 4. Using Bayes instead of JSE.")
      jse <- FALSE
    }
  }

  # get the prior means
  prior.means <- prior.means %||% getGlobalMeans(
    obj = original.x,
    targets = targets,
    assay = assay
  )

  is.atac_or_rna <- assay %in% c("atac", "rna")
  input.fun <- if (jse) {
    mean
  } else if (is.atac_or_rna) {
    atac_fun
  } else {
    median
  }

  input.assay <- if (is.atac_or_rna) {
    assays(original.x)$counts
  } else {
    flogit(assays(original.x)$Beta) # make sure we are with betas or we will double flogit
  }

  # bin the input
  bin.mat <- getBinMatrix(
    mat = cbind(input.assay, prior.means),
    genloc = rowRanges(x),
    chr = chr,
    res = res,
    FUN = input.fun,
    genome = genome
  )

  # shrink the bins using a James-Stein Estimator
  x.shrink <- t(apply(bin.mat$x, 1, .shrinkRow, jse, targets))

  # drop things that are zeroes as global means
  # this can and does crop up in resampling when you have something sparse
  # for instance single-cell data...
  # the correlation will break otherwise
  zeroes <- bin.mat$x[, "globalMean"] == 0
  if (any(zeroes)) {
    bin.mat$gr <- bin.mat$gr[!zeroes, ]
    x.shrink <- x.shrink[!zeroes, ]
    bin.mat$x <- bin.mat$x[!zeroes, ]
  }

  list(
    gr = bin.mat$gr,
    x = x.shrink[, colnames(x)],
    gmeans = bin.mat$x[, "globalMean"]
  )
}

.shrinkRow <- function(x, jse, targets) {
  x.samps <- x[!names(x) %in% "globalMean"]
  x.prior.m <- x["globalMean"]
  if (jse) {
    .jse(x = x.samps, grand.mean = x.prior.m, targets = targets)
  } else {
    .ebayes(x = x.samps, prior = x.prior.m, targets = targets)
  }
}

# helper function for summary when JSE == FALSE
atac_fun <- function(x) {
  sqrt(mean(x)) * length(x)
}

# helper functions for computing shrunken means
.ebayes <- function(x, prior = NULL, targets = NULL) {
  C <- if (is.null(targets)) sd(x) else sd(x[targets])
  prior + C * (x - prior)
}

#' James-Stein estimator
#' @param x input vector of binned means across samples
#' @param grand.mean The global mean across samples
#' @param targets Samples to shrink towards
#'
#' \eqn{\hat{\theta}_{JS+} = \left(1 - \frac{(m - 3)\sigma^2}{||\textbf{y} - \nu||^2}\right)}
#' @keywords internal
.jse <- function(x, grand.mean = NULL, targets = NULL) {
  m <- length(x)
  yv.norm <- sum((x - grand.mean)^2)
  sdsq <- if (is.null(targets)) sd(x)^2 else sd(x[targets])^2

  c <- 1 - (((m - 3) * sdsq) / yv.norm)
  grand.mean + c * (x - grand.mean)
}
