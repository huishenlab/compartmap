#' Plots A/B compartment estimates on a per chromosome basis
#'
#' Plot A/B compartments bins
#'
#' @param x      The matrix obejct returned from getCompartments
#' @param chr    Chromosome to subset to for plotting
#' @param what   Which metadata column to plot
#' @param main    Title for the plot
#' @param ylim      Y-axis limits (default is -1 to 1)
#' @param unitarize   Should the data be unitarized?
#' @param reverse    Reverse the sign of the PC values?
#' @param top.col    Top (pos. PC values) chromatin color to be plotted
#' @param bot.col    Bottom (neg. PC values) chromatin color to be plotted
#' @param with.ci    Whether to plot confidence intervals
#' @param filter    Whether to filter eigenvalues close to zero (default: TRUE)
#' @param filter.min.eigen  Minimum absolute eigenvalue to include in the plot
#' @param median.conf Plot the median confidence estimate across the chromosome?
#'
#' @import    GenomicRanges
#'
#' @export
#'
#' @return    A plot of inferred A/B compartments
#'
#' @examples
#'
#' library(GenomicRanges)
#'
#' # Generate random genomic intervals of 1-1000 bp on chr1-22
#' # Modified from https://www.biostars.org/p/225520/
#' random_genomic_int <- data.frame(chr = rep("chr14", 100))
#' random_genomic_int$start <- apply(random_genomic_int, 1, function(x) {
#'   round(runif(1, 0, getSeqLengths(chr = x)[[1]]), 0)
#' })
#' random_genomic_int$end <- random_genomic_int$start + runif(1, 1, 1000)
#' random_genomic_int$strand <- "*"
#'
#' # Generate random counts
#' counts <- rnbinom(1000, 1.2, 0.4)
#'
#' # Build random counts for 10 samples
#' count.mat <- matrix(sample(counts, nrow(random_genomic_int) * 10, replace = FALSE), ncol = 10)
#' colnames(count.mat) <- paste0("sample_", seq(1:10))
#'
#' # Bin counts
#' bin.counts <- getBinMatrix(count.mat, makeGRangesFromDataFrame(random_genomic_int), chr = "chr14", genome = "hg19")
#'
#' # Calculate correlations
#' bin.cor.counts <- getCorMatrix(bin.counts)
#'
#' # Get A/B signal
#' absignal <- getABSignal(bin.cor.counts)
#'
#' # Plot the A/B signal
#' par(mar = c(1, 1, 1, 1))
#' par(mfrow = c(1, 1))
#' plotAB(absignal, what = "pc")
plotAB <- function(x, chr = NULL, what = "score", main = "", ylim = c(-1, 1),
                   unitarize = FALSE, reverse = FALSE, top.col = "deeppink4",
                   bot.col = "grey50", with.ci = FALSE, filter = TRUE,
                   filter.min.eigen = 0.02, median.conf = FALSE) {
  # what are we plotting
  # what <- match.arg(what)
  # (no, don't do that)
  if (!what %in% names(mcols(x))) stop(what, " is not among names(mcols(x))")

  # check if plotting CI
  # FIXME: refactor this
  if (with.ci) {
    if (is(x, "GRanges")) {
      if (!is.null(chr)) x <- keepSeqlevels(x, chr, pruning.mode = "coarse")
      if (("conf.est" %in% names(mcols(x)))) {
        if (filter) {
          x <- x[abs(as(mcols(x)[what], "matrix")) > filter.min.eigen, ]
        }
        x.mat <- as(mcols(x)[what], "matrix")
        if (unitarize) x.mat <- .unitarize(x.mat)
        x.mat <- as.numeric(x.mat)
        if (reverse) x.mat <- -x.mat
        n <- length(x.mat)
        col <- rep(top.col, n)
        col[x.mat < 0] <- bot.col
        par(mar = c(1, 5, 1, 1))
        par(mfrow = c(2, 1))
        barplot(x.mat,
          ylim = ylim, bty = "n", xlab = "", ylab = "Eigenvector",
          border = col, col = col, main = main
        )

        barplot(x$conf.est, ylim = c(0, 1), ylab = "Compartment confidence estimate")
        if (median.conf) abline(h = median(x$conf.est), col = "red", lty = 2, lwd = 3)
      } else {
        message("conf.est isn't found in the mcols() of the input")
        stop("Run the compartmentCI() first.")
      }
    } else {
      stop("Input needs to be a GRanges object.")
    }
  } else {
    if (is(x, "GRanges")) {
      if (!is.null(chr)) x <- keepSeqlevels(x, chr, pruning.mode = "coarse")
      x <- as(mcols(x)[, what], "matrix")
    }
    if (unitarize) x <- .unitarize(x)
    if (filter) x <- x[abs(x) > filter.min.eigen]
    x <- as.numeric(x)
    if (reverse) x <- -x

    n <- length(x)
    col <- rep(top.col, n)
    col[x < 0] <- bot.col
    barplot(x,
      ylim = ylim, bty = "n", xlab = "", ylab = "Eigenvector",
      border = col, col = col, main = main
    )
  }
}


# helper fn
.unitarize <- function(x, medianCenter = TRUE) {
  if (medianCenter) x <- x - median(x, na.rm = TRUE)
  bad <- is.na(x)
  x[!bad] <- x[!bad] / sqrt(sum(x[!bad]^2))
  n.bad <- sum(bad)
  if (n.bad > 0) {
    message(
      sprintf("[.unitarize] %i missing values were ignored.\n", n.bad)
    )
  }
  x
}
