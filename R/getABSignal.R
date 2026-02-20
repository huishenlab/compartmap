#' Calculate Pearson correlations of smoothed eigenvectors
#' 
#' This function is used to generate a list x to be passed to getABSignal
#'
#' @param x      A list object from getCorMatrix
#' @param squeeze    Whether squeezing was used (implies Fisher's Z transformation)
#' @param assay What kind of assay are we working on ("array", "atac", "array")
#' @param genome The genome to use for gene-density-based sign correction
#' 
#' @return    A list x to pass to getABSignal
#' 
#' @import    SummarizedExperiment
#' 
#' @export 
#' 
#' @examples 
#' 
#' library(SummarizedExperiment)
#' library(BiocSingular)
#' 
#' #Generate random genomic intervals of 1-1000 bp on chr1-22
#' #Modified from https://www.biostars.org/p/225520/
#' random_genomic_int <- data.frame(chr = rep("chr14", 100))
#' random_genomic_int$start <- apply(random_genomic_int, 1, function(x) {
#'   round(runif(1, 0, getSeqLengths(getGenome("hg19"), chr = x)[[1]]), 0)
#' })

#' random_genomic_int$end <- random_genomic_int$start + runif(1, 1, 1000)
#' random_genomic_int$strand <- "*"
#' 
#' #Generate random counts
#' counts <- rnbinom(1000, 1.2, 0.4)
#' 
#' #Build random counts for 10 samples
#' count.mat <- matrix(sample(counts, nrow(random_genomic_int) * 10, replace = FALSE), ncol = 10)
#' colnames(count.mat) <- paste0("sample_", seq(1:10))
#' 
#' #Bin counts
#' bin.counts <- getBinMatrix(
#'   count.mat,
#'   makeGRangesFromDataFrame(random_genomic_int),
#'   chr = "chr14",
#'   genome = "hg19"
#' )
#'
#' #Calculate correlations
#' bin.cor.counts <- getCorMatrix(bin.counts)
#' 
#' #Get A/B signal
#' absignal <- getABSignal(bin.cor.counts)

getABSignal <- function(
  x,
  squeeze = FALSE,
  assay = c("rna", "atac", "array"),
  genome = c("hg19", "hg38", "mm9", "mm10")
) {
  assay <- match.arg(assay)
  gr <- x$gr

  message("Calculating eigenvectors.")
  pc <- getSVD(x$binmat.cor, sing.vec = "right")
  if (squeeze) pc <- ifisherZ(pc)

  message("Smoothing eigenvector.")
  gr$pc <- meanSmoother(pc)
  message("Done smoothing.")

  if (flipSign(gr, genome)) gr$pc <- -gr$pc
  gr$compartments <- extractOpenClosed(gr, assay = assay)
  return(gr)
}

flipSign <- function(gr, genome) {
  tx.gr <- getGenome(genome, "tx")
  gene_count <- countOverlaps(gr, tx.gr)
  open <- gr$pc > 0
  sum(gene_count[open]) < sum(gene_count[!open])
}

#' Get the open and closed compartment calls based on sign of singular values
#'
#' @param gr Input GRanges with associated mcols that represent singular values
#' @param cutoff Threshold to define open and closed states
#' @param assay The type of assay we are working with
#'
#' @return A vector of binary/categorical compartment states
#' @import SummarizedExperiment
#' @importFrom methods is
#' @export
#'
#' @examples
#'
#' dummy <- matrix(rnorm(10000), ncol = 25)
#' sing_vec <- getSVD(dummy, k = 1, sing.vec = "right")
#'
extractOpenClosed <- function(
  gr,
  cutoff = 0,
  assay = c("rna", "atac", "array")
) {
  # check for input to be GRanges
  if (!is(gr, "GRanges")) stop("Input needs to be a GRanges.")
  if (!("pc" %in% names(mcols(gr)))) stop("Need to have an mcols column be named 'pc'.")

  assay <- match.arg(assay)
  is.atac_or_rna <- assay %in% c("atac", "rna")
  is.open <- .isCompartmentOpen(is.atac_or_rna, gr$pc, cutoff)
  ifelse(is.open, "open", "closed")
}

# Check if a compartment is open based on assay type and eigenvalue
#
# For ATAC/RNA:
# eigen < cutoff - closed
# eigen > cutoff - open
#
# For methylation the logic is flipped:
# eigen < cutoff - open
# eigen > cutoff - closed
.isCompartmentOpen <- function(is.atac_or_rna, eigen, cutoff) {
  (is.atac_or_rna & eigen > cutoff) | (!is.atac_or_rna & eigen < cutoff)
}
