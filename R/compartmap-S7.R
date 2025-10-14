#' CompartmentCall class (experimental)
#'
#' An S7 class to hold a compartment call with its metadata and analysis
#' methods. This is primarily a parent class from which the more used
#' `CompartmapCall` and `MultiCompartmentCall` are derived.
#'
#' @param pc The singular values from a compartment call
#' @param res The binning resolution used
#' @param gr The GRanges of the bins
#' @param name An identifier for this copmartment call. Important to set if you
#' are making a MultiCompartmentCall
#' @param unitarized Whether the singular values have been unitarized
#'
#' @import S7
#' @importFrom data.table data.table
#' @export
CompartmentCall <- new_class(
  "CompartmentCall",
  properties = list(
    name = class_character,
    gr = methods::getClass("GRanges"),
    dt = methods::getClass("data.table"),
    res = class_numeric,
    unitarized = class_logical
  ),
  constructor = function(pc, res, gr, name = NULL, unitarized = FALSE) {
    dt <- data.table(pc = as.vector(pc))[, .(n = .I, pc, name)]
    new_object(
      S7_object(),
      name = name %||% shQuote(substitute(gr), "cmd2"),
      gr = granges(gr),
      dt = dt,
      res = res,
      unitarized = unitarized
    )
  }
)
S4_register(CompartmentCall)

#' Get singular values as a data.table
#' @export
DF <- new_generic("DF", "ccall")
method(DF, CompartmentCall) <- function(ccall) {
  ccall@dt[]
}

#' @name GRanges
#'
#' Get GRanges of the CompartmentCall
method(granges, CompartmentCall) <- function(x) {
  x@gr
}

#' Get the resolution of the CompartmentCall
#' @export
resolution <- new_generic("resolution", "ccall")
method(resolution, CompartmentCall) <- function(ccall) {
  ccall@res
}

#' Check if the CompartmentCall was unitarized
#' @export
is_unitarized <- new_generic("is_unitarized", "ccall")
method(is_unitarized, CompartmentCall) <- function(ccall) {
  ccall@unitarized
}

#' Get the name of the CompartmentCall object
#' @export
get_name <- new_generic("get_name", "ccall")
method(get_name, CompartmentCall) <- function(ccall) {
  ccall@name
}

#' Unitarize the singular values in a CompartmentCall or MultiCompartmentCall
#'
#' @export
unitarize <- new_generic("unitarize", "ccall")
method(unitarize, CompartmentCall) <- function(ccall, medianCenter = TRUE) {
  stopifnot("object is already unitarized" = isFALSE(ccall@unitarized))

  dt <- ccall@dt
  if ("name" %in% colnames(dt)) {
    ccall@dt <- ccall@dt[, .(n, pc = .unitarize(pc)), by = name][, .(n, pc, name)]
  } else {
    ccall@dt <- ccall@dt[, .(n, pc = .unitarize(pc))]
  }

  ccall@unitarized <- TRUE
  ccall
}

#' Flip the singular values signs in a CompartmentCall
#'
#' @export
flip <- new_generic("flip", "ccall")
method(flip, CompartmentCall) <- function(ccall) {
  ccall@dt <- ccall@dt[, .(n, pc = -pc, name)]
  ccall
}

#' Plot singular values of a CompartmentCall object
#'
#' @importFrom ggplot2 ggplot geom_line scale_y_continuous
#' @export
`plot.compartmap::CompartmentCall` <- function(x) {
  ggplot(x@dt, aes(x = n, y = pc)) +
    geom_line(linewidth = 1) +
    scale_y_continuous(limits = c(-1, 1))
}

