#' Denoising of Covariance matrix using Random Matrix Theory
#'
#' @name estRMT
#'
#' @details
#' This method takes in data as a matrix object. It then
#' fits a marchenko pastur density to eigenvalues of the correlation matrix. All
#' eigenvalues above the cutoff are retained and ones below the cutoff are
#' replaced such that the trace of the correlation matrix is 1 or non-significant
#' eigenvalues are deleted and diagonal of correlation matrix is changed to 1.
#' Finally, correlation matrix is converted to covariance matrix. This function
#' was taken and modified from the covmat package (https://github.com/cran/covmat)
#' which has since been deprecated on CRAN.
#'
#' @importFrom Matrix nearPD
#' @importFrom RMTstat dmp qmp
#'
#' @param  R input matrix
#' @param  Q ratio of rows/size. Can be supplied externally or fit using data
#' @param  cutoff takes two values max/each. If cutoff is max, Q is fitted and
#'          cutoff for eigenvalues is calculated. If cutoff is each, Q is set to
#'          row/size. Individual cutoff for each eigenvalue is calculated and used
#'          for filteration.
#' @param eigenTreat takes 2 values, average/delete. If average then the noisy
#'        eigenvalues are averged and each value is replaced by average. If delete
#'        then noisy eigenvalues are ignored and the diagonal entries of the
#'        correlation matrix are replaced with 1 to make the matrix psd.
#' @param numEig number of eigenvalues that are known for variance calculation.
#'        Default is set to 1. If numEig = 0 then variance is assumed to be 1.
#' @return A denoised RMT object
#'
#' @examples
#' rand_cor_mat <- cor(matrix(rnorm(100), nrow = 10))
#' denoised_rand_cor_mat <- estRMT(rand_cor_mat)$cov
#'
#' @author Rohit Arora
#'
#' @export
#'
estRMT <- function(R, Q = NA, cutoff = c("max", "each"),
                   eigenTreat = c("average", "delete"), numEig = 1) {
  .data <- as.matrix(R)
  T <- nrow(.data)
  M <- ncol(.data)
  if (T < M) stop("Does not work when nrow < ncol")

  if (!is.na(Q)) if (Q < 1) stop("Does not work for Q<1")
  cutoff <- cutoff[1]
  if (!cutoff %in% c("max", "each")) stop("Invalid cutoff")
  if (cutoff == "each") Q <- T / M

  eigenTreat <- eigenTreat[1]
  if (!eigenTreat %in% c("average", "delete")) stop("Invalid eigenTreat option")

  if (numEig < 0) stop("Number of eigenvalues must be non-negative")

  # eigenvalues can be negative. To avoid this e need a positive-definite matrix
  S <- cov(.data)
  S <- as.matrix(Matrix::nearPD(S)$mat)
  D <- diag(diag(S))
  C <- cov2cor(S)

  # Marchenko Pastur density is defined for eigenvalues of correlation matrix
  eigen.C <- eigen(C, symmetric = T)
  lambdas <- eigen.C$values
  sigma.sq <- mean(lambdas)

  sigma.sq <- 1 - sum(head(lambdas, numEig)) / M

  # minimize log-likelihood.
  loglik.marpas <- function(theta, sigma.sq) {
    Q <- theta
    val <- sapply(lambdas, function(x) RMTstat::dmp(x, svr = Q, var = sigma.sq))
    val <- val[val > 0]
    ifelse(is.infinite(-sum(log(val))), .Machine$double.xmax, -sum(log(val)))
  }

  if (is.na(Q) && cutoff != "each") {
    lb <- 1
    ub <- max(T / M, 5)
    starts <- seq(lb, ub, length.out = 50)
    # this would be a logical place to use BiocParallel::bpapply
    fit.marpas <- do.call(
      rbind,
      lapply(
        starts,
        function(start) {
          optim(
            par = start,
            fn = loglik.marpas,
            method = "L-BFGS-B",
            lower = lb,
            upper = ub,
            sigma.sq = sigma.sq
          )
        }
      )
    )
    idx <- grep("CONVERGENCE", unlist(fit.marpas[, "message"]))
    vals <- fit.marpas[idx, c("par", "value")] # wtf is going on here
    Q <- unlist(vals[which.min(vals[, "value"]), "par"])
  }

  lambda.max <- RMTstat::qmp(1, svr = Q, var = sigma.sq)
  # now that we have a fit. lets denoise eigenvalues below the cutoff

  if (cutoff == "max") {
    idx <- which(lambdas > lambda.max)
  } else if (cutoff == "each") {
    cutoff.each <- sapply(2:length(lambdas), function(i) {
      eigr <- lambdas[i:M]
      mean(eigr) * (1 + (M - i + 1) / T + 2 * sqrt((M - i + 1) / T))
    })
    idx <- c(1, 1 + which(lambdas[-1] > cutoff.each))
  }

  if (length(idx) == 0) {
    return(S)
  }

  val <- eigen.C$values[idx]
  vec <- eigen.C$vectors[, idx, drop = FALSE]
  sum <- 0
  for (i in 1:ncol(vec)) sum <- sum + val[i] * vec[, i] %*% t(vec[, i])

  # trace of correlation matrix is 1. Use this to determine all the remaining
  # eigenvalues

  lambdas.cleaned <- c()
  clean.C <- if (eigenTreat == "average") {
    lambdas.cleaned <- c(val, rep(1, M))
    sum + sum(eigen.C$values[-idx]) / M * diag(rep(1, M))
  } else if (eigenTreat == "delete") {
    lambdas.cleaned <- c(val, rep(0, M))
    diag(sum) <- 1
    sum
  }

  # convert correlation to covariance matrix and return
  clean.S <- D^0.5 %*% clean.C %*% D^0.5
  fit <- list(
    cov = clean.S, Q = Q, var = sigma.sq, eigVals = lambdas,
    eigVals.cleaned = lambdas.cleaned, lambdascutoff = lambda.max
  )

  class(fit) <- "RMT"
  fit
}

#' Wrapper to denoise a correlation matrix using a Random Matrix Theory approach
#'
#' @name getDenoisedMatrix
#'
#' @param obj SummarizedExperiment object with rowRanges for each feature and colnames
#' @param res The resolution desired (default is a megabase 1e6)
#' @param chr Which chromosome to perform the denoising
#' @param genome Which genome (default is hg19)
#' @param iter How many iterations to perform denoising
#' @param targets Samples/cells to shrink towards
#' @param prior.means The means of the bin-level prior distribution (default will compute them for you)
#' @param assay What assay type this is ("rna", "atac")
#'
#' @return A denoised correlation matrix object for plotting with plotCorMatrix
#'
#' @import SummarizedExperiment
#' @import GenomicRanges
#' @import scales
#'
#' @export
#'
#' @examples
#' data("k562_scrna_chr14", package = "compartmap")
#' denoised_cor_mat <- getDenoisedCorMatrix(k562_scrna_chr14, genome = "hg19", assay = "rna")
getDenoisedCorMatrix <- function(obj, res = 1e6, chr = "chr14",
                                 genome = c("hg19", "hg38", "mm9", "mm10"),
                                 iter = 2, targets = NULL, prior.means = NULL,
                                 assay = c("rna", "atac", "array")) {
  ## this is a wrapper to give back a denoised correlation matrix to plot
  # match the assay args
  assay <- match.arg(assay)

  # match the genome if given
  genome <- match.arg(genome)

  # double check the obj class is compatible
  if (!checkAssayType(obj)) stop("Input needs to be a SummarizedExperiment")

  # subset to selected chromosome(s)
  obj <- keepSeqlevels(obj, chr, pruning.mode = "coarse")

  # get the prior means
  if (is.null(prior.means)) {
    prior.means <- getGlobalMeans(
      obj = obj, targets = targets,
      assay = assay
    )
  } else {
    message("Assuming the prior means passed were derived from the full sample set.")
  }

  # shrink bins
  message("Shrinking bins with the JSE.")
  bin.mat <- shrinkBins(
    x = obj, original.x = obj, prior.means = prior.means,
    chr = chr, res = res, targets = targets,
    jse = TRUE, assay = assay, genome = genome
  )
  # get the raw correlation matrix
  cor.mat <- getCorMatrix(binmat = bin.mat, squeeze = FALSE)
  # denoise with RMT
  message("Denoising the correlation matrix using RMT.")
  cor.mat.denoise <- compartmap::estRMT(cor.mat$binmat.cor)$cov
  # iterate?
  if (iter >= 2) {
    for (i in 2:iter) {
      message("Iterative denoising. Iteration: ", i)
      cor.mat.denoise <- compartmap::estRMT(cor.mat.denoise)$cov
    }
  }

  # rescale
  cor.mat.denoise <- scales::rescale(cor.mat.denoise,
    to = c(0, 1)
  )
  # tidy up
  colnames(cor.mat.denoise) <- rownames(cor.mat.denoise) <- as.character(granges(cor.mat$gr.cor))

  return(cor.mat.denoise)
}

#' Plot a denoised correlation matrix
#'
#' @name plotCorMatrix
#'
#' @param denoised.cor.mat The denoised correlation matrix object from getDenoisedMatrix
#' @param midpoint The midpoint for the coloring (default is 0.3)
#' @param return.plot.obj Whether to return the ggplot object
#' @param uppertri Whether to keep the upper triangle of the matrix
#' @param lowertri Whether to keep the lower triangle of the matrix
#'
#' @return Either a ggplot object or plot
#'
#' @import ggplot2
#' @import reshape2
#'
#' @export
#'
#' @examples
#' dummy <- matrix(rnorm(10000), ncol = 25)
#' set.seed(1000)
#' my_plot <- plotCorMatrix(dummy, return.plot.obj = TRUE)
plotCorMatrix <- function(denoised.cor.mat,
                          midpoint = 0.3,
                          return.plot.obj = FALSE,
                          uppertri = FALSE,
                          lowertri = FALSE) {
  ## upper tri
  if (uppertri) {
    denoised.cor.mat[upper.tri(denoised.cor.mat)] <- 0
  }
  if (lowertri) {
    denoised.cor.mat[lower.tri(denoised.cor.mat)] <- 0
  }
  diag(denoised.cor.mat) <- 1
  ## melt for plotting
  cor.mat.melt <- reshape2::melt(denoised.cor.mat)
  ## plot
  p <- ggplot(cor.mat.melt, aes(x = Var2, y = Var1, fill = value)) +
    geom_raster() +
    scale_fill_gradient2(
      low = "white", mid = "white", high = "red3", midpoint = midpoint,
      name = "Correlation"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank(),
      plot.margin = grid::unit(c(0, 0, 0, 0), "in")
    )

  if (return.plot.obj) {
    return(p)
  } else {
    p
  }
}
