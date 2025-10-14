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

#' CompartmapCall class (experimental)
#'
#' An S7 class to hold a single-sample or grouped compartment call with its
#' metadata and analysis methods. This can take the GRanges output of
#' `scCompartments(group = TRUE)`. For multiple single-cell level compartment
#' inferences use MultiCompartmentCall.
#'
#' @param gr The GRanges output of scCompartments or getArrayCompartments
#' containing the 'pc' column
#' @param res The binning resolution used
#' @param name An identifier for this compartment call. Important to set if you
#' are making a MultiCompartmentCall.
#' @param unitarized Whether the singular values have been unitarized
#'
#' @export
CompartmapCall <- new_class(
  "CompartmapCall",
  parent = CompartmentCall,
  constructor = function(gr, res, name = NULL, unitarized = FALSE) {
    dt <- data.table(pc = gr$pc)[, `:=`(n = .I)][]
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
S4_register(CompartmapCall)

#' MultiCompartmentCall class
#'
#' An S7 class to hold multiple compartment calls at the same resolution with
#' their metadata and analysis methods. This can take a set of single
#' `CompartmentCall` or `CompartmapCall` objects.
#'
#' @param ccalls A list of CompartmapCalls to combine
#' @param name An identifier for this set of compartment calls
#' @param unitarized Whether the singular values have been unitarized
#' @param unitarize Whether to unitarize the singular values for each of the inputs calls
#'
#' @export
#' @importFrom data.table rbindlist dcast
#' @export
MultiCompartmentCall <- new_class(
  "MultiCompartmentCall",
  parent = CompartmentCall,
  properties = list(
    colnames = class_character,
    mat = new_S3_class(c("matrix", "array"))
  ),
  constructor = function(ccalls, name, unitarized = FALSE, unitarize = FALSE) {
    all_same <- function(prop_get, err) {
      unique_property <- unique(unlist(lapply(ccalls, prop_get)))
      if (length(unique_property) != 1) stop(err)
      unique_property
    }

    unique_res <- all_same(resolution, "All resolutions must be the same")
    unique_gr <- all_same(granges, "All calls must be either unitarized or non-unitarized")

    if (unitarize) {
      if (unitarized) {
        message("All singular values already unitarized")
      } else {
        ccalls <- lapply(seq_along(ccalls), function(i) {
          if (!is_unitarized(ccalls[[i]])) {
            return(unitarize(ccalls[[i]]))
          }
          ccalls[[i]]
        })
      }
      unitarized <- TRUE
    }

    all_unitarized <- all_same(is_unitarized, "All calls must be either unitarized or non-unitarized")

    dt <- rbindlist(lapply(ccalls, function(i) {
      DF(i)[, name := get_name(i)][]
    }))
    mat <- as.matrix(dcast(dt, n ~ name, value.var = "pc")[, -1])

    new_object(
      S7_object(),
      name = name,
      gr = unique_gr[[1]],
      dt = dt,
      res = unique_res,
      unitarized = unitarized,
      colnames = unlist(lapply(ccalls, get_name)),
      mat = mat
    )
  },
  validator = function(self) {
    if (length(self@colnames) != length(unique(self@colnames))) {
      "All names must be unique"
    }
  }
)
S4_register(MultiCompartmentCall)

#' Compute agreement between compartment calls
#'
#' Calculates the proportion of calls with the same sign for every pair of
#' calls in a MultiCompartmentCall object.
#'
#' @export
agr <- new_generic("agr", "mcall")
method(agr, MultiCompartmentCall) <- function(mcall) {
  agreement(mcall@mat)
}

#' Compute correlation between compartment calls
#'
#' Calculates Pearson correlation for every pair of calls in a
#' MultiCompartmentCall object.
#'
#' @importFrom stats cor
#' @export
corr <- new_generic("corr", "x")
method(corr, MultiCompartmentCall) <- function(x) {
  cor(x@mat)
}

#' Plot singular values from a MultiCompartmentCall object
#'
#' @importFrom ggplot2 ggplot geom_line scale_y_continuous
#' @export
`plot.compartmap::MultiCompartmentCall` <- function(x) {
  ggplot(x@dt, aes(x = n, y = pc, color = name)) +
    geom_line(linewidth = 1) +
    scale_y_continuous(limits = c(-1, 1))
}
