#' Generate bins for A/B compartment estimation
#'
#' Generate bins across a user defined chromosome for A/B compartment estimation.
#' A/B compartment estimation can be used for non-supported genomes if chr.end is set.
#'
#' This function is used to generate a list object to be passed to getCorMatrix
#'
#' @param mat       A p x n matrix where p (rows) = loci and n (columns) = samples/cells
#' @param genloc    GRanges object that contains corresponding genomic locations of the loci
#' @param chr       Chromosome to be analyzed
#' @param chr.start Starting position (in bp) to be analyzed
#' @param chr.end   End position (in bp) to be analyzed
#' @param res       Binning resolution (in bp)
#' @param FUN       Function to be used to summarize information within a bin
#' @param genome    Genome corresponding to the input data ("hg19", "hg38", "mm9", "mm10")
#'
#' @return    A list object to pass to getCorMatrix
#' @import    SummarizedExperiment
#' @export
#'
#' @examples
#' library(GenomicRanges)
#'
#' # Generate random genomic intervals of 1-1000 bp on chr1-22
#' # Modified from https://www.biostars.org/p/225520/
#' genome.gr <- getGenome("hg19")
#' random_genomic_int <- data.frame(chr = rep("chr14", 100))
#' random_genomic_int$start <- apply(random_genomic_int, 1, function(x) {
#'   round(runif(1, 0, getSeqLengths(genome.gr, chr = x)[[1]]), 0)
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
getBinMatrix <- function(
  mat,
  genloc,
  chr = "chr1",
  chr.start = 0,
  chr.end = NULL,
  res = 100000,
  FUN=sum,
  genome = c("hg19", "hg38", "mm9", "mm10")
) {

  if (any(is.na(mat))){
    stop("Matrix must not contain NAs")
  }
  if (nrow(mat) != length(genloc)){
    stop("Provided GRanges must have length equal to the matrix number of rows")
  }

  chr.end <- chr.end %||% getSeqLengths(getGenome(genome), chr = chr)
  gr.bin <- .makeBins(chr.start, chr.end, res, chr)
  ids <- findOverlaps(genloc, gr.bin, select = "first")

  binCount <- length(gr.bin)
  flog.debug("%d bins created...", binCount)

  mat.bin <- apply(mat, 2, function(x) {
    .summarizeBins(x, binCount, ids, FUN)
  })
  colnames(mat.bin) <- colnames(mat)

  # Subset the non-empty bins
  wh <- rowSums(mat.bin) != 0
  mat.bin <- mat.bin[wh,]
  gr.bin  <- gr.bin[wh]

  if (nrow(mat.bin) < 2) stop("There are not enough non-empty bins to continue...")

  list(gr = gr.bin, x = mat.bin)
}

.makeBins <- function(start, end, res, chr) {
  starts <- seq(start, end, by = res)
  ends <- c(starts[-1], end) - 1L
  GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start = starts, end = ends)
  )
}

.summarizeBins <- function(matCol, binCount, ids, fun) {
  zvec <- rep(0, binCount)
  a <- tapply(matCol, INDEX = ids, FUN = fun)  # Summarize with user defined function
  zvec[as.numeric(names(a))] <- a
  zvec
}
