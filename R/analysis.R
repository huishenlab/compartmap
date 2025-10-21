#' Compute sign agreement of a matrix of singular values
#' @param mat A matrix of compartment calls where each column is a single-cell
#' or grouped set of singular values
#' @importFrom utils combn
#' @export
agreement <- function(mat) {
  cmbs <- combn(1:ncol(mat), 2)
  ncols <- ncol(mat)
  agrmat <- matrix(nrow = ncols, ncol = ncols)
  for (cidx in seq(ncol(cmbs))) {
    i <- cmbs[, cidx][1]
    j <- cmbs[, cidx][2]
    res <- .agreement(mat[, i], mat[, j])
    agrmat[i, j] <- res
    agrmat[j, i] <- res
  }
  dimnames(agrmat) <- list(colnames(mat), colnames(mat))
  diag(agrmat) <- 1
  agrmat
}

.agreement <- function(v1, v2) {
  b1 <- v1 > 0
  b2 <- v2 > 0
  sum(b1 == b2) / length(v1)
}
