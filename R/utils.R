#' Check if the assay is a SummarizedExperiment
#'
#' @param obj Input object
#' @return NULL
#' @importFrom methods is
#' @keywords internal
#'
#' @examples
#' data("k562_scrna_chr14", package = "compartmap")
#' compartmap:::verifySE(k562_scrna_chr14)
verifySE <- function(obj) {
  # helper function to check the class of an object
  if (!is(obj, "SummarizedExperiment")) {
    stop("Input needs to be a SummarizedExperiment")
  }
}


#' Throw error if assay does not contain coordinates
#'
#' @param obj Input object
#'
#' @return NULL
#' @keywords internal
#'
#' @examples
#' data("k562_scrna_chr14", package = "compartmap")
#' compartmap:::verifyCoords(k562_scrna_chr14)
verifyCoords <- function(obj) {
  # helper function to check the class of an object
  if (length(seqinfo(rowRanges(obj))) == 0) {
    stop(paste(
      "The SummarizedExperiment you have provided has no coordinates.\n",
      "Compartment extraction will fail.\n",
      "Please provide rowRanges with genomic coordinates for the object."
    ))
  }
}

#' Check that the input SummarizedExperiment object has the right assays
#'
#' @param se Input SummarizedExperiment object
#' @param assay The assay type
#'
#' @return Error if the right assay type is not present, NULL if it is
#' @keywords internal
verifyAssayNames <- function(se, assay) {
  reqName <- switch(assay,
    rna = "counts",
    atac = "counts",
    array = "Beta",
    bisulfite = "Beta",
    stop(shQuote(assay), " is unsupported")
  )
  if (!reqName %in% assayNames(se)) {
    stop("The 'assays' slot should contain ", shQuote(reqName), " for ", assay, " data")
  }
}

#' Get the assay names from a SummarizedExperiment object
#'
#' @param se Input SummarizedExperiment object
#'
#' @return The names of the assays
#' @export
#'
#' @examples
#' data("k562_scrna_chr14", package = "compartmap")
#' getAssayNames(k562_scrna_chr14)
getAssayNames <- function(se) {
  # helper function to check the assay slot names
  names(assays(se))
}

#' Helper function: squeezed logit
#'
#' @param p       a vector of values between 0 and 1 inclusive
#' @param sqz     the amount by which to 'squeeze', default is .000001
#'
#' @return        a vector of values between -Inf and +Inf
#'
#' @examples
#'
#' p <- runif(n = 1000)
#' summary(p)
#'
#' sqz <- 1 / (10**6)
#' x <- flogit(p, sqz = sqz)
#' summary(x)
#'
#' all(abs(p - fexpit(x, sqz = sqz)) < sqz)
#' all(abs(p - fexpit(flogit(p, sqz = sqz), sqz = sqz)) < sqz)
#'
#' @export
flogit <- function(p, sqz = 0.000001) {
  midpt <- 0.5
  deflate <- 1 - (sqz * midpt)
  if (any(p > 1 | p < 0, na.rm = TRUE)) stop("Values of p outside (0,1) detected.")
  squoze <- ((p - midpt) * deflate) + midpt
  return(log(squoze / (1 - squoze)))
}

#' Helper function: expanded expit
#'
#' @param x       a vector of values between -Inf and +Inf
#' @param sqz     the amount by which we 'squoze', default is .000001
#'
#' @return        a vector of values between 0 and 1 inclusive
#'
#' @examples
#'
#' x <- rnorm(n = 1000)
#' summary(x)
#'
#' sqz <- 1 / (10**6)
#' p <- fexpit(x, sqz = sqz)
#' summary(p)
#'
#' all((abs(x - flogit(p)) / x) < sqz)
#' all(abs(x - flogit(fexpit(x))) < sqz)
#'
#' @export
fexpit <- function(x, sqz = 0.000001) {
  midpt <- .5
  squoze <- exp(x) / (1 + exp(x))
  inflate <- 1 / (1 - (sqz * midpt))
  ((squoze - midpt) * inflate) + midpt
}

#' Get the chromosomes from an object
#'
#' @param obj Input SummarizedExperiment object
#'
#' @return A character vector of chromosomes present in an object
#' @import SummarizedExperiment
#'
#' @examples
#' data("k562_scrna_chr14", package = "compartmap")
#' getChrs(k562_scrna_chr14)
#' @export
getChrs <- function(obj) {
  # get the chromosomes present in the object
  return(unique(as.character(seqnames(obj))))
}

#' Remove bootstrap estimates that failed
#'
#' This can happen if the correlation between the bins and eigenvector fails
#' theoretically we can recover these but need an additional utility to find
#' consensus
#'
#' @param obj Input list object with elements 'pc' and 'gr'
#' @return A filtered list object
#' @export
#' @keywords internal
removeEmptyBoots <- function(obj) {
  Filter(Negate(anyNA), obj)
}

#' Get a GRanges object from bundled compartmap genomes
#'
#' @param genome The desired genome to use ("hg19", "hg38", "mm9", "mm10")
#' @param type The type of data - full genome or open sea regions
#'
#' @return Granges of the genome
#'
#' @examples
#' hg19 <- getGenome(genome = "hg19")
#'
#' @export
#' @keywords internal
getGenome <- function(
  genome = c("hg19", "hg38", "mm9", "mm10"),
  type = "genome"
) {
  genome.name <- match.arg(genome) |> tryCatch(error = function(e) {
    e <- gsub("'arg'", "'genome'", e)
    msg <- paste0(e, "Only human and mouse genomes are supported for the time being.")
    stop(msg)
  })
  gr <- switch(type,
    genome = paste0(genome.name, ".gr"),
    openseas = paste0("openSeas.", genome.name)
  )
  return(get(gr))
}

#' Get the seqlengths of a chromosome from a given genome's GRanges
#'
#' The goal for this function is to eliminate the need to lug around
#' large packages when we only want seqlengths for things.
#'
#' @param genome.gr A GRanges object of the genome (from `getGenome()`)
#' @param chr What chromosome to extract the seqlengths of
#'
#' @return The seqlengths of a specific chromosome
#'
#' @importFrom GenomeInfoDb seqlengths seqlevels
#' @import GenomicRanges
#'
#' @examples
#' hg19.chr14.seqlengths <- getSeqLengths(getGenome('hg19'), chr = "chr14")
#'
#' @export
getSeqLengths <- function(genome.gr, chr = "chr14") {
  sl <- seqlengths(genome.gr)[chr]
  if (is.na(sl)) {
    genome.build <- unique(genome(genome.gr))
    msg <- paste(
      chr, "not found in seqlevels of", genome.build,
      "- check that the 'genome' and 'chr' arguments are correct"
    )
    stop(msg)
  }
  return(sl)
}

#' Get chunked sets of row-wise or column-wise indices of a matrix
#'
#' @name getMatrixBlocks
#'
#' @param mat Input matrix
#' @param chunk.size The size of the chunks to use for coercion
#' @param by.row Whether to chunk in a row-wise fashion
#' @param by.col Whether to chunk in a column-wise fashion
#'
#' @return A set of chunked indices
#'
#' @examples
#' # make a sparse binary matrix
#' library(Matrix)
#' m <- 100
#' n <- 1000
#' mat <- round(matrix(runif(m * n), m, n))
#' mat.sparse <- Matrix(mat, sparse = TRUE)
#'
#' # get row-wise chunks of 10
#' chunks <- getMatrixBlocks(mat.sparse, chunk.size = 10)
#'
#' @export
getMatrixBlocks <- function(
  mat,
  chunk.size = 1e5,
  by.row = TRUE,
  by.col = FALSE
) {
  message("Using chunk size: ", chunk.size)
  if (by.row) {
    message("Breaking into row chunks.")
    return(split(1:nrow(mat), ceiling(seq_along(1:nrow(mat)) / chunk.size)))
  }

  # assumes column-wise chunking
  message("Breaking into column chunks.")
  return(split(1:ncol(mat), ceiling(seq_along(1:ncol(mat)) / chunk.size)))
}

#' Convert a sparse matrix to a dense matrix in a block-wise fashion
#'
#' @name sparseToDenseMatrix
#'
#' @param mat Input sparse matrix
#' @param blockwise Whether to do the coercion in a block-wise manner
#' @param by.row Whether to chunk in a row-wise fashion
#' @param by.col Whether to chunk in a column-wise fashion
#' @param chunk.size The size of the chunks to use for coercion
#' @param parallel Whether to perform the coercion in parallel
#' @param cores The number of cores to use in the parallel coercion
#'
#' @return A dense matrix of the same dimensions as the input
#'
#' @import Matrix
#' @importFrom parallel mclapply
#' @importFrom methods as
#'
#'
#' @examples
#' # make a sparse binary matrix
#' library(Matrix)
#' m <- 100
#' n <- 1000
#' mat <- round(matrix(runif(m * n), m, n))
#' mat.sparse <- Matrix(mat, sparse = TRUE)
#'
#' # coerce back
#' mat.dense <- sparseToDenseMatrix(mat.sparse, chunk.size = 10)
#'
#' # make sure they are the same dimensions
#' dim(mat) == dim(mat.dense)
#'
#' # make sure they are the same numerically
#' all(mat == mat.dense)
#'
#' @export
sparseToDenseMatrix <- function(
  mat,
  blockwise = TRUE,
  by.row = TRUE,
  by.col = FALSE,
  chunk.size = 1e5,
  parallel = FALSE,
  cores = 2
) {
  if (isFALSE(blockwise)) {
    return(as(mat, "matrix"))
  }

  # do block-wise reconstruction of matrix
  chunks <- getMatrixBlocks(mat,
    chunk.size = chunk.size,
    by.row = by.row, by.col = by.col
  )

  if (by.row & parallel) {
    return(do.call("rbind", mclapply(chunks, function(r) {
      return(as(mat[r, ], "matrix"))
    }, mc.cores = cores)))
  }

  if (by.row & !parallel) {
    return(do.call("rbind", lapply(chunks, function(r) {
      return(as(mat[r, ], "matrix"))
    })))
  }

  # assumes column-wise conversion
  if (by.col & parallel) {
    return(do.call("cbind", mclapply(chunks, function(r) {
      return(as(mat[, r], "matrix"))
    }, mc.cores = cores)))
  }

  return(do.call("cbind", lapply(chunks, function(r) {
    return(as(mat[, r], "matrix"))
  })))
}

#' Import and optionally summarize a bigwig at a given resolution
#'
#' @name importBigWig
#'
#' @param bw Path a bigwig file
#' @param bins Optional set of bins as a GRanges to summarize the bigwig to
#' @param summarize Whether to perform mean summarization
#' @param genome Which genome is the bigwig from ("hg19", "hg38", "mm9", "mm10")
#'
#' @return SummarizedExperiment object with rowRanges corresponding to summarized features
#'
#' @import SummarizedExperiment
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqlengths seqlevels seqlevelsStyle<-  keepSeqlevels keepStandardChromosomes
#'
#' @export

importBigWig <- function(
  bw,
  bins = NULL,
  summarize = FALSE,
  genome = c("hg19", "hg38", "mm9", "mm10")
) {
  # read in the bigwig
  bw.raw <- rtracklayer::import(bw)
  # coerce to UCSC style seqlevels
  seqlevelsStyle(bw.raw) <- "UCSC"
  if (!is.null(bins)) {
    seqlevelsStyle(bins) <- "UCSC"
  }
  # it is now a GRanges object
  if (any(is.na(seqlengths(bw.raw)))) stop("Imported bigwig does not have seqlengths")
  ## only supporting human and mouse for now
  if (genome %in% c("hg19", "hg38")) {
    species <- "Homo_sapiens"
  } else {
    species <- "Mus_musculus"
  }
  bw.sub <- keepStandardChromosomes(bw.raw, species = species, pruning.mode = "coarse")
  if (!is.null(bins)) {
    bins <- keepSeqlevels(bins, value = seqlevels(bw.sub), pruning.mode = "coarse")
  }
  if (summarize) {
    # make sure it's sorted
    bw.sub <- sort(bw.sub)
    # this assumes seqlengths exist...
    # this also assumes some bins exist
    if (is.null(bins)) stop("Specify bins as GRanges with tileGenome")
    bw.score <- GenomicRanges::coverage(bw.sub, weight = "score")
    bw.bin <- GenomicRanges::binnedAverage(bins, bw.score, "ave_score")
    # cast to a SummarizedExperiment to bin them
    bw.se <- SummarizedExperiment(
      assays = S4Vectors::SimpleList(counts = as.matrix(mcols(bw.bin)$ave_score)),
      rowRanges = granges(bw.bin)
    )
    colnames(bw.se) <- as.character(bw)
    return(bw.se)
  }
  return(bw.sub)
}

#' Generate function to filter rows/columns with NAs exceeding a threshold
#'
#' @details
#' Since removing NAs from rows vs columns only differs by whether rowMeans or
#' colMeans is used, and by where the comma goes in the subset operation,
#' code repetition can be avoided by consolidating these operations.
#' This `cleanAssay` function can generate two functions to remove NA's from
#' rows and columns using the `by` argument based on which it selects the
#' appropriate 'mean' and subset functions. This maintains the clarity of
#' having the operation in the function name when used: `cleanAssayRows` and
#' `cleanAssayCols`.
#' @param by Whether to filter by rows or columns
#'
#' @return A function to filter assay rows/columns
#' @keywords internal
cleanAssay <- function(by = c("row", "col")) {
  by <- match.arg(by)
  if (by == "row") {
    mean.fun <- rowMeans
    subset.fun <- function(se, toKeep) {
      se[toKeep, ]
    }
  } else {
    mean.fun <- colMeans
    subset.fun <- function(se, toKeep) {
      se[, toKeep]
    }
  }

  function(se, na.max = 0.8, assay = c("array", "bisulfite")) {
    assay <- match.arg(assay)
    assay.data <- switch(assay,
      array = assays(se)$Beta,
      bisulfite = assays(se)$counts
    )
    toKeep <- mean.fun(is.na(assay.data)) < na.max
    subset.fun(se, toKeep)
  }
}

#' Remove rows with NAs exceeding a threshold. See `cleanAssay()`
#'
#' @param se Input SummarizedExperiment object
#' @param na.max The maximum number of NAs allowed as a fraction
#' @param assay The type of assay we are working with
#'
#' @return A filtered matrix
#'
#' @examples
#' if (requireNamespace("minfi", quietly = TRUE)) {
#'   data("array_data_chr14", package = "compartmap")
#'   compartmap:::cleanAssayRows(array.data.chr14, assay = "array")
#' }
#' @keywords internal
cleanAssayRows <- cleanAssay(by = "row")

#' Remove columns/cells/samples with NAs exceeding a threshold. See `cleanAssay()`
#'
#' @param se Input SummarizedExperiment object
#' @param na.max The maximum number of NAs allowed as a fraction
#' @param assay The type of assay we are working with
#'
#' @return A filtered matrix
#'
#' @examples
#' if (requireNamespace("minfi", quietly = TRUE)) {
#'   data("array_data_chr14", package = "compartmap")
#'   compartmap:::cleanAssayCols(array.data.chr14, assay = "array")
#' }
#' @keywords internal
cleanAssayCols <- cleanAssay(by = "col")

#' Filter to open sea CpG loci
#'
#' @name filterOpenSea
#'
#' @param obj Input SummarizedExperiment or GRanges object
#' @param genome Which genome to filter
#' @param other GRanges of open sea regions (TODO)
#'
#' @return Filtered to open sea CpG loci
#' @import SummarizedExperiment
#' @importFrom methods is
#' @importFrom utils data
#' @export
#'
#' @examples
#' if (requireNamespace("minfi", quietly = TRUE)) {
#'   data("array_data_chr14", package = "compartmap")
#'   opensea <- filterOpenSea(array.data.chr14, genome = "hg19")
#' }
#'
#' @export
filterOpenSea <- function(
  obj,
  genome = c("hg19", "hg38", "mm10", "mm9"),
  other = NULL
) {
  stopifnot("'obj' needs to be a GRanges or SummarizedExperiment" = is(obj, "GRanges") | is(obj, "SummarizedExperiment"))

  # get the desired open sea loci given the genome GRanges
  openseas.genome <- other %||% getGenome(genome, type = "openseas")
  stopifnot("The 'other' input needs to be a GRanges of open sea regions" = is(openseas.genome, "GRanges"))

  # Subset by overlaps
  message("Filtering to open sea CpG loci...")
  # subset to just CpG loci if CpH or rs probes still exist
  obj <- obj[grep("cg", rownames(obj)), ]
  subsetByOverlaps(obj, openseas.genome)
}

#' Gather open sea CpG from a GRanges of CpG islands
#'
#' @description This function accepts a GRanges input of CpG islands that can
#' be derived from UCSC table browser and rtracklayer::import(yourbed.bed)
#'
#' @name filterOpenSea
#'
#' @param gr Input GRanges of CpG islands
#'
#' @return GRanges object that can be used with filterOpenSea()
#' @import rtracklayer
#' @import GenomicRanges
#' @export
#'
#' @examples
#' #cpgi <- rtracklayer::import(system.file("inst/extdata/mm10_cpgi.bed", package = "compartmap"))
#' #opensea_cpg <- getOpenSeas(cpgi)
getOpenSeas <- function(gr) {
  resorts <- trim(resize(gr, width(gr) + 8000, fix = "center"))
  openSeas <- subset(gaps(resorts), strand == "*")
  return(openSeas)
}
