#' Transform/normalize counts using TF-IDF
#'
#' @details
#' This function and its helpers were modeled after or taken from:
#' - http://andrewjohnhill.com/images/posts/2019-5-6-dimensionality-reduction-for-scatac-data/analysis.html
#' - https://divingintogeneticsandgenomics.rbind.io/post/clustering-scatacseq-data-the-tf-idf-way/
#'
#' @param mat n x p input matrix (n = samples/cells; p = rna counts)
#' @param scale.factor Scaling factor for the term-frequency (TF)
#' @param cap The maximum expression count used for TF-IDF. Default of 1
#' binarizes the matrix. A `cap` value greater than 1 will cap counts at that
#' value.
#'
#' @return A TF-IDF transformed matrix of the same dimensions as the input
#'
#' @import Matrix
#'
#' @examples
#' m <- 1000
#' n <- 100
#' mat <- round(matrix(runif(m * n), m, n))
#' # Input needs to be a tall matrix
#' tfidf <- transformTFIDF(mat)
#'
#' @export
transformTFIDF <- function(mat, scale.factor = 1e5, cap = 1) {
  if (!is(mat, "matrix") & !is(mat, "Matrix")) {
    stop("Input needs to be a matrix.")
  }

  # this assumes n x p matrix (e.g. a wide matrix)
  # check and transpose as needed
  # input matrix is tall
  mat.capped <- Matrix(mat, sparse = TRUE)
  if (dim(mat)[1] < dim(mat)[2]) {
    mat.capped <- t(Matrix(mat, sparse = TRUE))
  }

  # binarize the matrix
  expr.limit <- ifelse(cap == 1, 0, cap)
  mat.capped@x[mat.capped@x > expr.limit] <- cap
  tf <- t(t(mat.capped) / Matrix::colSums(mat.capped))           # compute term-frequency
  tf@x <- log1p(tf@x * scale.factor)                             # scale
  idf <- log(1 + ncol(mat.capped) / Matrix::rowSums(mat.capped)) # inverse-document frequency smooth
  tfidf <- .tfidf(tf, idf)                                       # transform

  # cast back to a matrix since things like UMAP don't like sparse matrices
  as.matrix(t(tfidf))
}

# helper function for TF-IDF transform
# modeled after http://andrewjohnhill.com/images/posts/2019-5-6-dimensionality-reduction-for-scatac-data/analysis.html
.tfidf <- function(tf, idf) {
  tf <- t(tf)
  tf@x <- tf@x * rep.int(idf, diff(tf@p))
  tf <- t(tf)
  return(tf)
}

#' Transform/normalize compartment calls using TF-IDF on HDF5-backed objects
#'
#' @name hdf5TFIDF
#'
#' @param h5 SummarizedExperiment object, DelayedMatrix, or a normal matrix
#' @param scale.factor Scaling factor for the term-frequency (TF)
#' @param return.dense Whether to return a dense, in memory matrix
#' @param return.se Whether to return the TF-IDF matrix as a new assay in the SummarizedExperiment
#'
#' @return A TF-IDF transformed matrix of the same dimensions as the input
#' @import Matrix
#' @import DelayedMatrixStats
#' @import DelayedArray
#' @import HDF5Array
#' @importFrom methods as is
#'
#' @examples
#'
#' m <- 1000
#' n <- 100
#' mat <- round(matrix(runif(m * n), m, n))
#' # Input needs to be a tall matrix
#' tfidf <- hdf5TFIDF(mat)
#'
#' @export
hdf5TFIDF <- function(h5, scale.factor = 1e5,
                      return.dense = FALSE,
                      return.se = FALSE) {
  # binarze
  if (is(h5, "SummarizedExperiment")) {
    if (!is(assay(h5), "DelayedMatrix")) {
      assay(h5) <- as(assay(h5), "HDF5Matrix") # coerce to hdf5 backing to work with any SE
    }
    assay(h5)[assay(h5) > 0] <- 1
    h5.mat <- assay(h5) # tall matrix
  }

  if (is(h5, "DelayedMatrix")) {
    # make the matrix tall if needed
    if (dim(h5)[1] < dim(h5)[2]) h5 <- t(h5)
    h5[h5 > 0] <- 1
    h5.mat <- h5
  }

  if (is(h5, "matrix")) {
    h5[h5 > 0] <- 1
    h5.mat <- as(h5, "HDF5Matrix")
  }

  message("Computing term frequency.")
  tf <- t(t(h5.mat) / DelayedMatrixStats::colSums2(h5.mat)) # term frequency
  tf <- log1p(tf * scale.factor)                            # scale

  message("Computing inverse document frequency.")
  idf <- log(1 + ncol(h5.mat) / DelayedMatrixStats::rowSums2(h5.mat)) # inverse document frequency

  # TODO: fix this ugliness...
  tf.mat <- as.matrix(tf)
  tf.sparse <- Matrix(tf.mat, sparse = TRUE) # cast the tf matrix back to a sparse matrix
  tf.sparse <- t(tf.sparse)                  # transpose for TF-IDF

  message("TF-IDF")
  tf.sparse@x <- tf.sparse@x * rep.int(idf, diff(tf.sparse@p)) # TF-IDF applied
  tf.sparse <- t(tf.sparse)                                    # transpose again

  # coerce back to dense matrix
  if (return.dense) {
    message("WARNING: This might blow up!")
    message("If you get a cholmod error: problem too large, set return.dense to FALSE.")
    message("You will get a sparse matrix returned instead.")
    return(as.matrix(tf.sparse))
  }
  if (return.se) {
    if (!is(h5, "SummarizedExperiment")) {
      return(tf.sparse)
    }
    message("Returning the TF-IDF matrix into the SummarizedExperiment.")
    assays(h5)$tfidf <- tf.sparse
    return(h5)
  }
  return(tf.sparse)
}
