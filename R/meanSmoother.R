#' Windowed mean smoother
#'
#' @name meanSmoother
#'
#' @param mat       Input data matrix: samples are columns, regions/loci are rows
#' @param k         Number of windows to use (default = 1 to smooth with 1
#'                  window on either side of a position)
#' @param iters     Number of iterations to smooth (default is 2)
#' @param delta     Convergence threshhold (overrides iter if > 0; default is 0)
#' @param weights   Weights, if using any (NULL)
#'
#' @importFrom stats median weighted.mean
#'
#' @return      Smoothed data matrix
#'
#' @examples
#' dummy <- matrix(rnorm(10000), ncol = 25)
#' smooth.dummy <- meanSmoother(dummy)
#' smooth.dummy <- meanSmoother(dummy, iters = 3)
#' smooth.dummy <- meanSmoother(dummy, delta = 1e-3)
#'
#' @export
meanSmoother <- function(mat, k = 1, iters = 2, delta = 0, weights = NULL) {
  if (k == 0) {
    message("Returning unsmoothed mat as 'k' = 0")
    return(mat)
  }

  stopifnot(length(mat) >= k)

  pos <- 0
  eps <- delta + 1
  weights <- weights %||% rep(1, length(mat))

  while (pos < iters & eps > delta) {
    mat0 <- mat
    pos <- pos + 1
    mat <- .meanSmoother.internal(mat0, weights = weights, k = k)
    eps <- median(abs(mat - mat0), na.rm = TRUE)
  }

  mat
}


# helper fn
# finds the first and last eligible positions to smooth as well as excess bins
# beyond eligible positions.
.meanSmoother.internal <- function(mat, weights, k) {
  n <- length(mat)
  y <- rep(NA, n)

  first_pos <- k + 1
  last_pos <- n - k
  excess_bins <- seq((last_pos + 1), n)

  # note that na.rm can create issues
  for (pos in first_pos:last_pos) {
    y[pos] <- .window.mean(mat, weights, pos, k = k)
  }

  # it is possible for windows to be 0 in these loops
  for (pos in 1:k) {
    y[pos] <- .window.mean(mat, weights, pos, k = pos - 1)
  }

  for (pos in excess_bins) {
    y[pos] <- .window.mean(mat, weights, pos, k = n - pos)
  }

  y
}


# helper fn
.window.mean <- function(mat, weights, pos, k) {
  if (k < 1) {
    return(mat[pos])
  }
  endpos <- pos + k
  startpos <- pos - k - 1
  stride <- seq(startpos, endpos)
  weighted.mean(mat[stride], weights = weights[stride], na.rm = TRUE)
}
