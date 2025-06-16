#' Plots A/B compartment estimates on a per chromosome basis
#'
#' Plot A/B compartments bins
#'
#' @param grAB              The GRanges object returned from scCompartments and getArrayABsignal
#' @param chr               Chromosome to subset to for plotting
#' @param main              Title for the plot
#' @param ylim              Y-axis limits (default is -1 to 1)
#' @param unitarize         Should the data be unitarized?
#' @param reverse           Reverse the sign of the PC values?
#' @param top.col           Top (pos. PC values) chromatin color to be plotted
#' @param bot.col           Bottom (neg. PC values) chromatin color to be plotted
#' @param with.ci           Whether to plot confidence intervals
#' @param filter            Whether to filter eigenvalues close to zero (default: TRUE)
#' @param filter.min.eigen  Minimum absolute eigenvalue to include in the plot
#' @param median.conf       Plot the median confidence estimate across the chromosome?
#'
#' @import GenomicRanges
#' @importFrom methods as is
#' @importFrom stats median
#' @importFrom graphics abline barplot par
#'
#' @return    A plot of inferred A/B compartments
#'
#' @export
#'
#' @examples
#' library(GenomicRanges)
#'
#' # Generate random genomic intervals of 1-1000 bp on chr1-22
#' # Modified from https://www.biostars.org/p/225520/
#' random_genomic_int <- data.frame(chr = rep("chr14", 100))
#' random_genomic_int$start <- apply(random_genomic_int, 1, function(x) {
#'   round(runif(1, 0, getSeqLengths(getGenome("hg19"), chr = x)[[1]]), 0)
#' })
#' random_genomic_int$end <- random_genomic_int$start + runif(1, 1, 1000)
#' random_genomic_int$strand <- "*"
#'
#' # Generate random counts
#' counts <- rnbinom(1000, 1.2, 0.4)
#'
#' # Build random counts for 10 samples
#' count.mat <- matrix(
#'   sample(counts, nrow(random_genomic_int) * 10, replace = FALSE),
#'   ncol = 10
#' )
#' colnames(count.mat) <- paste0("sample_", seq(1:10))
#'
#' # Bin counts
#' bin.counts <- getBinMatrix(
#'   count.mat,
#'   makeGRangesFromDataFrame(random_genomic_int),
#'   chr = "chr14",
#'   genome = "hg19"
#' )
#'
#' # Calculate correlations
#' bin.cor.counts <- getCorMatrix(bin.counts)
#'
#' # Get A/B signal
#' absignal <- getABSignal(bin.cor.counts)
#'
#' # Plot the A/B signal
#' plotAB(absignal)
plotAB <- function(
  grAB,
  chr = NULL,
  main = "",
  ylim = c(-1, 1),
  unitarize = FALSE,
  reverse = FALSE,
  top.col = "deeppink4",
  bot.col = "grey50",
  with.ci = FALSE,
  filter = TRUE,
  filter.min.eigen = 0.02,
  median.conf = FALSE
) {
  stopifnot("'grAB' is not a GRanges object" = is(grAB, "GenomicRanges"))
  mcolnames <- names(mcols(grAB))

  if (with.ci && !("conf.est" %in% mcolnames)) {
    stop("conf.est isn't found in the mcols() of the input - run the compartmentCI() first.")
  }

  if (!is.null(chr)) grAB <- keepSeqlevels(grAB, chr, pruning.mode = "coarse")
  mat.AB <- as(mcols(grAB)[, "pc"], "matrix")
  if (unitarize) mat.AB <- .unitarize(mat.AB)
  if (filter) mat.AB <- mat.AB[abs(mat.AB) > filter.min.eigen]
  if (reverse) mat.AB <- -mat.AB
  mat.AB <- as.numeric(mat.AB)

  n <- length(mat.AB)
  col <- rep(top.col, n)
  col[mat.AB < 0] <- bot.col
  if (with.ci) {
    par(mar = c(1, 5, 1, 1), mfrow = c(2, 1))
    .barplotAB(mat.AB, ylim, col, main)
    barplot(grAB$conf.est, ylim = c(0, 1), ylab = "Compartment confidence estimate")
    if (median.conf) abline(h = median(grAB$conf.est), col = "red", lty = 2, lwd = 3)
  } else {
    .barplotAB(mat.AB, ylim, col, main)
  }
}

# helper fn
.unitarize <- function(x, medianCenter = TRUE) {
  if (medianCenter) x <- x - median(x, na.rm = TRUE)
  NAs <- is.na(x)
  x[!NAs] <- x[!NAs] / sqrt(sum(x[!NAs]^2))
  na.count <- sum(NAs)
  if (na.count > 0) message(sprintf("[.unitarize] %i missing values were ignored.\n", na.count))
  x
}

.barplotAB <- function(mat.AB, ylim, col, main) {
  barplot(
    mat.AB,
    ylim = ylim,
    bty = "n",
    xlab = "",
    ylab = "Eigenvector",
    border = col,
    col = col,
    main = main
  )
}
