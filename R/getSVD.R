#' Compute the SVD of a matrix using irlba
#'
#' @name getSVD
#'
#' @param matrix A p x n input matrix
#' @param k Number of singular vectors to return
#' @param sing.vec Whether to return the right or left singular vector
#'
#' @return A singular vector or matrix with sign corresponding to positive values
#'
#' @importFrom stats cor
#' @importFrom BiocSingular IrlbaParam runSVD
#' @export
#'
#' @examples
#'
#' dummy <- matrix(rnorm(10000), ncol=25)
#' sing_vec <- getSVD(dummy, k = 1, sing.vec = "right")
getSVD <- function(matrix, k = 1, sing.vec = c("left", "right")) {
  sing.vec <- match.arg(sing.vec)

  matrix <- .centerMatrix(matrix)
  V <- .getUV(matrix, k, sing.vec)
  V * sqrt(length(V)) # Chromosome length normalization
}

# Centre the matrix subtracting the row means from each row
#' @keywords internal
.centerMatrix <- function(mat) {
  t(scale(t(mat), center = TRUE, scale = FALSE))
}

#' Run SVD and get left or right singular vectors
#'
#' @param matrix A p x n input matrix
#' @param k Number of singular vectors to return
#' @param sing.vec Whether to return the right or left singular vector
#'
#' @return left or right singular vectors
#'
#' @keywords internal
.getUV <- function(matrix, k, sing.vec) {
  SVD <- runSVD(matrix, k = k, BSPARAM = IrlbaParam())
  switch(sing.vec, left = SVD$u, right = SVD$v)
}
