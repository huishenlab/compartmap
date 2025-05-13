#' Fisher's Z transformation
#'
#' `fisherZ` returns (squeezed) Fisher's Z transformed Pearson's r
#'
#' @param cormat Pearson correlation matrix
#'
#' @return Fisher Z transformed Pearson correlations
#' @export
#'
#' @examples
#'
#' #Generate a random binary (-1, 1) matrix
#' (mat <- matrix(sample(c(1,-1), 25, replace = TRUE), ncol = 5))
#'
#' #Correct matrix diag
#' diag(mat) <- 1
#'
#' #Transform
#' fisherZ(mat)
fisherZ <- function(cormat) {
  if (any(cormat == 1)) cormat <- .squeezeit(cormat)
  atanh(cormat)
}

# Helper function to squeeze binary matrix for transformation
.squeezeit <- function(cormat) {
  cormat * 0.999999
}

#' Inverse Fisher's Z transformation
#'
#' `ifisherZ` returns the inverse (squeezed) Fisher's Z transformed
#' Pearson's r.
#'
#' @param zmat matrix of Fisher's Z transformed Pearson correlations or an
#' eignevector
#'
#' @return Back transformed Fisher's Z Pearson correlations
#' @export
#'
#' @examples
#'
#' # Generate a random binary (-1, 1) matrix
#' (mat <- matrix(sample(c(1,-1), 25, replace = TRUE), ncol = 5))
#'
#' # Correct matrix diag
#' diag(mat) <- 1
#'
#' # Transform
#' (mat.transform <- fisherZ(mat))
#'
#' #Back transform
#' ifisherZ(mat.transform)
ifisherZ <- function(zmat) {
  tanh(zmat)
}
