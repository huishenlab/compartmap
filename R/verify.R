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
  reqName <- switch(
    assay,
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
