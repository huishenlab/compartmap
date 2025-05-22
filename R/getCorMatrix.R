#' Calculate Pearson correlations of a binned matrix
#'
#' This function is used to generate a list object to be passed to getABSignal
#'
#' @param binmat      A binned matrix list object from getBinMatrix
#' @param squeeze    Whether to squeeze the matrix for Fisher's Z transformation
#'
#' @return    A list object to pass to getABSignal
#'
#' @import    GenomicRanges
#' @importFrom stats cor
#'
#' @export
#'
#' @examples
#'
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
#' count.mat <- matrix(sample(counts, nrow(random_genomic_int) * 10, replace = FALSE), ncol = 10)
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
#' getCorMatrix(bin.counts)
getCorMatrix <- function(binmat, squeeze = FALSE) {
  message("Calculating correlations...")
  # bind back up the global means and shrunken bins
  binmat$x <- cbind(binmat$x, binmat$gmeans)
  binmat.cor <- suppressWarnings(cor(t(binmat$x)))

  if (squeeze) {
    binmat.cor <- fisherZ(binmat.cor)
  }
  message("Done...")
  list(gr.cor = binmat$gr, binmat.cor = binmat.cor)
}
