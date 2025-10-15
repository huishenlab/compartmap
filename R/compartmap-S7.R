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
    dt <- data.table(pc = as.vector(pc))[, .(n = .I, pc, name = name)]
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

#' @export
method(`[`, CompartmentCall) <- function(x, i = NULL) {
  i <- i %||% seq_len(nrow(x@mat))
  x@dt <- x@dt[i]
  x@dt[, n := .I][]
  x@gr <- x@gr[i]
  x
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
#' @param x CompartmentCall or CompartmapCall object
#' @param label_coords Label the x-axis with genomic coordinates. Uses a
#' numeric index when set to `FALSE`. Using coordinate labels can severely
#' crowd the x-axis, especially with Kb-resolution calls.
#' @param res The resolution to round the genomic coordinates to (kilobase:
#' "kb" or megabase: "mb")
#'
#' @importFrom ggplot2 ggplot geom_line scale_y_continuous theme element_text
#' @export
`plot.compartmap::CompartmentCall` <- function(x, label_coords = FALSE, res = "kb", linewidth = 0.5) {
  if (label_coords) {
    p <- ggplot(x@dt[, .(n, pc, coord = grscale(x@gr, res))], aes(x = coord, y = pc, group = 1)) +
      geom_line(linewidth = linewidth) +
      scale_y_continuous(limits = c(-1, 1)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    return(p)
  }
  ggplot(x@dt, aes(x = n, y = pc)) +
    geom_line(linewidth = linewidth) +
    scale_y_continuous(limits = c(-1, 1))
}

#' Print CompartmentCall
method(print, CompartmentCall) <- function(x) {
  message(.print_CompartmentCall(x))
}

.print_CompartmentCall <- function(x) {
  properties <- props(x)
  class_type <- gsub("compartmap::", "", class(x)[1])

  sprintf(
    "<%s object>
  @name        : %s
  @res         : %s
  @gr          : GRanges with %d bins
  @dt          : data.table of compartment calls (n = bin index, pc = singular values)
  @unitarized  : %s",
    class_type,
    x@name,
    .resolution(x@res),
    length(x@gr),
    x@unitarized
  )
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
    dt <- data.table(pc = gr$pc)[, `:=`(n = .I, name = name)][, .(n, pc, name)]
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
    unique_gr <- all_same(granges, "All GRanges must contain the same ranges")

    all_unitarized <- unique(unlist(lapply(ccalls, is_unitarized)))
    if (unitarized & !all_unitarized) {
      stop("Not all calls are unitarized - unitarize all inputs or run with `unitarize = TRUE`")
    }

    if (unitarize) {
      if (all_unitarized) {
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

.resolution <- function(res) {
  if (res == 1e5) {
    paste(res / 1e5, "Kb")
  } else if (res == 1e6) {
    paste(res / 1e6, "Mb")
  } else {
    stop("Unsupported resolution")
  }
}

method(print, MultiCompartmentCall) <- function(x) {
  msg <- message(
    .print_CompartmentCall(x),
    sprintf("\n  @mat         : %d bins x %d samples", nrow(x@mat), ncol(x@mat))
  )
}

method(`[`, MultiCompartmentCall) <- function(x, i = NULL) {
  i <- i %||% seq_len(nrow(x@mat))
  x@gr <- x@gr[i]
  x@mat <- x@mat[i, ]
  x@dt <- x@dt[n %in% i]
  x@dt[, n := seq_len(.N), by = name]
  x
}

#' @export
method(`[`, MultiCompartmentCall) <- function(x, i = NULL, j = NULL) {
  i <- i %||% seq_len(nrow(x@mat))
  j <- j %||% seq_len(ncol(x@mat))

  if (length(i) == 1 & length(j) == 1) {
    stop("Subsetting to a single value is not allowed")
  }
  stopifnot("i exceeds the row count of x" = i <= nrow(x@mat))

  if (is.numeric(j)) {
    if (length(j) > ncol(x@mat)) {
      stop("j exceeds the column count of x")
    }
    subset_names <- x@dt[, unique(name)][j]
  } else if (is.character(j)) {
    stopifnot("One or more columns not found in the data" = all(j %in% colnames(x@mat)))
    subset_names <- j
  } else {
    stop("Unsupported type for subsetting columns")
  }

  x@gr <- x@gr[i]
  x@mat <- x@mat[i, j, drop = FALSE]
  x@dt <- x@dt[n %in% i & name %in% subset_names]
  x
}

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
#' @param x MultiCompartmentCall
#' @param label_coords Label the x-axis with genomic coordinates. Uses a
#' numeric index when set to `FALSE`. Using coordinate labels can severely
#' crowd the x-axis, especially with Kb-resolution calls.
#' @param res The resolution to round the genomic coordinates to (kilobase:
#' "kb" or megabase: "mb")
#'
#' @importFrom ggplot2 ggplot geom_line scale_y_continuous theme element_text
#' @export
`plot.compartmap::MultiCompartmentCall` <- function(x, label_coords = FALSE, res = "mb", linewidth = 0.5) {
  if (label_coords) {
    p <- ggplot(
      x@dt[, .(n, pc, name, coord = grscale(x@gr, res))],
      aes(x = coord, y = pc, color = name, group = name)
    ) +
      geom_line(linewidth = linewidth) +
      scale_y_continuous(limits = c(-1, 1)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    return(p)
  }
  ggplot(x@dt, aes(x = n, y = pc, color = name)) +
    geom_line(linewidth = linewidth) +
    scale_y_continuous(limits = c(-1, 1))
}

grscale <- function(gr, res) {
  scale_factor <- switch(tolower(res), kb = list(1e5, "Kb"), mb = list(1e6, "Mb"))
  start_scaled <- start(gr) / as.numeric(scale_factor[1])
  end_scaled <- round(end(gr) / as.numeric(scale_factor[1], 4))
  paste0(seqlevels(gr), ":", start_scaled, "-", end_scaled, " ", scale_factor[2])
}
