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
#' @importFrom data.table data.table .I .N :=
#' @export
CompartmentCall <- new_class(
  "CompartmentCall",
  properties = list(
    name = class_character,
    gr = methods::getClass("GRanges"),
    dt = methods::getClass("data.table"),
    res = class_numeric,
    unitarized = class_logical,
    seqinfo = methods::getClass("Seqinfo")
  ),
  constructor = function(pc, res, gr, name = NULL, unitarized = FALSE) {
    dt <- data.table(pc = as.vector(pc))[, .(n = .I, pc, name = name)]
    new_object(
      S7_object(),
      name = name %||% shQuote(substitute(gr), "cmd2"),
      gr = granges(gr),
      dt = dt,
      res = res,
      unitarized = unitarized,
      seqinfo = selectMethod('seqinfo', "GRanges")(gr)
    )
  }
)
S4_register(CompartmentCall)

#' Get the @dt slot from a CompartmentCall object
#'
#' @param x A CompartmentCall object
#'
#' @export
DF <- new_generic("DF", "x", function(x) {
  S7_dispatch()
})
method(DF, CompartmentCall) <- function(x) {
  x@dt[]
}

#' Subset rows of a CompartmentCall object
#'
#' @param x A CompartmentCall object
#' @param i Row indices to subset
#'
#' @export
`[.compartmap::CompartmentCall` <- function(x, i = NULL) {
  i <- i %||% seq_len(nrow(x@mat))
  x@dt <- x@dt[i]
  x@dt[, n := .I][]
  x@gr <- x@gr[i]
  x
}

#' Get GRanges of the CompartmentCall
#'
#' @param x A CompartmentCall object
#'
#' @export
method(granges, CompartmentCall) <- function(x) {
  x@gr
}

#' Get Seqinfo of the CompartmentCall
#'
#' @param x A CompartmentCall object
#'
#' @export
method(seqinfo, CompartmentCall) <- function(x) {
  x@seqinfo
}


#' Get 'seqlevels' of the CompartmentCall
#'
#' @param x A CompartmentCall object
#'
#' @export
method(seqlevels, CompartmentCall) <- function(x) {
  selectMethod('seqlevels', 'GRanges')(x)
}

#' Subset the CompartmentCall object by chromosome
#'
#' @param x A CompartmentCall object
#' @param chr A string vector of chromosomes to subset to
#'
#' @export
subset_chr <- new_generic("subset_chr", "x", function(x, chr) {
  S7_dispatch()
})
method(subset_chr, CompartmentCall) <- function(x, chr) {
  ind <- as.vector(seqnames(x@gr) %gin% chr)
  x[which(ind)]
}

#' Find overlaps between CompartmentCall objects
#'
#' @param subject A CompartmentCall object
#' @param query A string vector of chromosomes to subset to
#'
#' @export
method(findOverlaps, list(CompartmentCall, CompartmentCall)) <- function(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = FALSE
) {
  query <- query@gr
  subject <- subject@gr
  callGeneric()
}

#' Get the resolution of the CompartmentCall
#'
#' @param x A CompartmentCall object
#'
#' @export
resolution <- new_generic("resolution", "x", function(x) {
  S7_dispatch()
})
method(resolution, CompartmentCall) <- function(x) {
  x@res
}

#' Check if the CompartmentCall was unitarized
#'
#' @param x A CompartmentCall object
#'
#' @export
is_unitarized <- new_generic("is_unitarized", "x", function(x) {
  S7_dispatch()
})
method(is_unitarized, CompartmentCall) <- function(x) {
  x@unitarized
}

#' Get the name of the CompartmentCall object
#'
#' @param x A CompartmentCall object
#'
#' @export
get_name <- new_generic("get_name", "x", function(x) {
  S7_dispatch()
})
method(get_name, CompartmentCall) <- function(x) {
  x@name
}

#' Unitarize the singular values in a CompartmentCall or MultiCompartmentCall
#'
#' @param x A CompartmentCall object
#' @param medianCenter Whether to center the singular values on their median
#'
#' @export
unitarize <- new_generic("unitarize", "x", function(x, medianCenter = TRUE) {
  S7_dispatch()
})
method(unitarize, CompartmentCall) <- function(x, medianCenter = TRUE) {
  stopifnot("object is already unitarized" = isFALSE(x@unitarized))

  dt <- x@dt
  if (inherits(x, "compartmap::MultiCompartmentCall")) {
    x@mat <- apply(x@mat, 2, .unitarize)
    x@dt <- x@dt[, .(n, pc = .unitarize(pc)), by = name][, .(n, pc, name)]
  } else {
    x@dt <- x@dt[, .(n, pc = .unitarize(pc), name)]
  }

  x@unitarized <- TRUE
  x
}

#' Flip the singular values signs in a CompartmentCall
#'
#' @param x A CompartmentCall object
#'
#' @export
flip <- new_generic("flip", "x", function(x) {
  S7_dispatch()
})
method(flip, CompartmentCall) <- function(x) {
  x@dt <- x@dt[, .(n, pc = -pc, name)]
  x
}

#' Fill missing genomic bins in CompartmentCalls using a reference GRanges
#'
#' Compartmap may drop genomic bins with insufficient data and the resulting
#' GRanges object may not have all the bins of the region it was run which
#' means using the same region and resolution on different inputs does not
#' guarantee the same output bins. Having different set of bins between
#' calls despite the same input regions and resolution prevents the creation of
#' MultiCompartmentCall objcets that expect the same GRanges across all input
#' CompartmentCall objects. `fill_missing()` adds the missing bins according to
#' a larger reference set of bins, filling missing data with NA. All
#' CompartmentCall objects bins must be present in the reference bins.
#'
#' @param x A CompartmentCall object
#'
#' @export
fill_missing <- new_generic("flip", "x", function(x, ref.gr) {
  S7_dispatch()
})
method(fill_missing, CompartmentCall) <- function(x, ref.gr) {
  ref_length <- length(ref.gr)
  stopifnot("Reference GRanges is not bigger than CompartmentCall object" = ref_length > length(x@gr))
  stopifnot("All CompartmentCall bins must be present in the reference GRanges" = all(x@gr %gin% ref.gr))

  ref_idx <- seq_len(ref_length)
  dt <- data.table(n = ifelse(ref.gr %gin% x@gr, ref_idx, NA))
  dt[!is.na(n), `:=`(pc = x@dt$pc, name = x@name)]
  dt[, n := .I][]
  x@gr <- ref.gr
  x@dt <- dt
  x
}

#' Plot singular values of a CompartmentCall object
#'
#' @param x CompartmentCall or CompartmapCall object
#' @param ... Placeholder for the `plot` generic - arguments have not effect
#' @param type Whether to plot the singular values as `"line"` or `"bar"` plots.
#' @param label_coords Label the x-axis with genomic coordinates. Uses a
#' numeric index when set to `FALSE`. Using coordinate labels can severely
#' crowd the x-axis, especially with Kb-resolution calls.
#' @param res The resolution to round the genomic coordinates to (kilobase:
#' "kb" or megabase: "mb")
#' @param width The width of the `geom_line` if `type = "line"` or the width
#' of the bar if `type = "bar"` in the plot
#' @param ylim Upper and lower bound for the y-axis
#'
#' @importFrom ggplot2 ggplot geom_line geom_bar scale_y_continuous theme element_text
#' @export
`plot.compartmap::CompartmentCall` <- function(
  x,
  ...,
  type = "line",
  label_coords = FALSE,
  res = "kb",
  width = 0.5,
  ylim = c(-1, 1)
) {
  . <- NULL
  pc <- NULL
  if (label_coords) {
    coord_pd <- x@dt[, .(n, pc, coord = grscale(x@gr, res))]
    p <- ggplot(coord_pd, aes(x = coord, y = pc, group = 1)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  } else {
    p <- ggplot(x@dt, aes(x = n, y = pc))
  }
  p <- p + scale_y_continuous(limits = ylim)

  switch(
    type,
    line = p + geom_line(linewidth = width),
    bar = p + geom_bar(stat = "identity", width = width)
  )
}

#' Print CompartmentCall
method(print, CompartmentCall) <- function(x, ...) {
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
      unitarized = unitarized,
      seqinfo = selectMethod("seqinfo", "GRanges")(gr)
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

    all_unitarized <- all(unique(unlist(lapply(ccalls, is_unitarized))))
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
    colorder <- sapply(ccalls, get_name)
    dt <- rbindlist(lapply(ccalls, function(i) {
      DF(i)[, name := get_name(i)][, name := factor(name, levels = colorder)][]
    }))
    mat <- as.matrix(dcast(dt, n ~ name, value.var = "pc")[, -1])

    new_object(
      S7_object(),
      name = name,
      gr = unique_gr[[1]],
      dt = dt,
      res = unique_res,
      unitarized = unitarized,
      seqinfo = selectMethod('seqinfo', "GRanges")(unique_gr[[1]]),
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
  fct <- res / 1e5
  if ((any(fct < 1) | any(fct >= 100))) stop("Unsupported resolution")
  ifelse(fct < 10, paste(fct, "Kb"), paste(fct / 10, "Mb"))
}

#' Print a CompartmentCall object
#'
#' @param x A CompartmentCall object
#'
#' @export
method(print, MultiCompartmentCall) <- function(x, ...) {
  column_label <- ifelse(inherits(x, "compartmap::SingleCellCompartmentCall"), "cells", "samples")
  msg <- message(
    .print_CompartmentCall(x),
    sprintf("\n  @mat         : %d bins x %d %s", nrow(x@mat), ncol(x@mat), column_label)
  )
}

#' Subset MultiCompartmentCall object
#'
#' @param x A CompartmentCall object
#' @param i Rows indices to subset
#' @param j Column indices or names to subset
#'
#' @export
`[.compartmap::MultiCompartmentCall` <- function(x, i = NULL, j = NULL) {
  i <- i %||% seq_len(nrow(x@mat))
  j <- j %||% seq_len(ncol(x@mat))

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
  x@dt[, n := seq_len(.N), by = name][]
  x
}


#' Compute agreement between compartment calls
#'
#' Calculates the proportion of calls with the same sign for every pair of
#' calls in a MultiCompartmentCall object.
#'
#' @param x A MultiCompartmentCall object
#'
#' @export
agr <- new_generic("agr", "x", function(x) {
  S7_dispatch()
})
method(agr, MultiCompartmentCall) <- function(x) {
  agreement(x@mat)
}

#' Compute correlation between compartment calls
#'
#' Calculates Pearson correlation for every pair of calls in a
#' MultiCompartmentCall object.
#'
#' @param x A MultiCompartmentCall object
#' @param ... Additional arguments to pass to `stats::cor()`
#'
#' @importFrom stats cor
#' @export
corr <- new_generic("corr", "x", function(x, ...) {
  S7_dispatch()
})
method(corr, MultiCompartmentCall) <- function(x, ...) {
  cor(x@mat, ...)
}

#' Plot singular values from a MultiCompartmentCall object
#'
#' @param x MultiCompartmentCall
#' @param ... Placeholder for the `plot` generic - arguments have not effect
#' @param type Whether to plot the singular values as `"line"` or `"bar"`
#' plots. Bar plots will be facted by the CompartmapCall object name while the
#' line plots are overlayed.
#' @param label_coords Label the x-axis with genomic coordinates. Uses a
#' numeric index when set to `FALSE`. Using coordinate labels can severely
#' crowd the x-axis, especially with Kb-resolution calls.
#' @param res The resolution to round the genomic coordinates to (kilobase:
#' "kb" or megabase: "mb")
#' @param width The width of the `geom_line` if `type = "line"` or the width
#' of the bar if `type = "bar"` in the plot
#' @param ylim Upper and lower bound for the y-axis
#'
#' @importFrom ggplot2 ggplot geom_line scale_y_continuous theme element_text
#' @export
`plot.compartmap::MultiCompartmentCall` <- function(
  x,
  ...,
  type = "line",
  label_coords = FALSE,
  res = "mb",
  width = 0.5,
  ylim = c(-1, 1)
) {
  pd <- x@dt
  x_axis <- "n"
  if (label_coords) {
    pd <- x@dt[, .(n, pc, name, coord = grscale(x@gr, res))]
    x_axis <- "coord"
  }

  p <- switch(
    type,
    line = {
      ggplot(pd, aes(x = .data[[x_axis]], y = pc, color = name, group = name)) +
        geom_line(linewidth = width)
    },
    bar = {
      ggplot(pd, aes(x = .data[[x_axis]], y = pc, group = name)) +
        geom_col(width = width) +
        facet_grid(rows = vars(name))
    }
  )

  if (label_coords) {
    p <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }
  p + scale_y_continuous(limits = ylim)
}

grscale <- function(gr, res) {
  scale_factor <- switch(tolower(res), kb = list(1e5, "Kb"), mb = list(1e6, "Mb"))
  start_scaled <- start(gr) / as.numeric(scale_factor[1])
  end_scaled <- round(end(gr) / as.numeric(scale_factor[1], 4))
  paste0(seqlevels(gr), ":", start_scaled, "-", end_scaled, " ", scale_factor[2])
}

#' SingleCellCompartmentCall class
#'
#' An S7 class to hold multiple single-cell level compartment calls at the same
#' resolution with their metadata and analysis methods. This takes a
#' RaggedExperiment from `scCompartments()` as its input.
#'
#' @param ccalls A RageedExperiment of single-cell compartment calls
#' @param name An identifier for this set of compartment calls
#' @param unitarized Whether the singular values have been unitarized
#' @param unitarize Whether to unitarize the singular values for each of the inputs calls
#'
#' @export
#' @importFrom data.table melt as.data.table
#' @export
SingleCellCompartmentCall <- new_class(
  "SingleCellCompartmentCall",
  parent = MultiCompartmentCall,
  properties = list(
    colnames = class_character,
    mat = new_S3_class(c("matrix", "array"))
  ),
  constructor = function(ccalls, res, name, unitarized = FALSE, unitarize = FALSE) {
    grlist <- condenseSE(ccalls)
    pcs <- lapply(grlist, function(i) {
      mcols(i)[, 'pc']
    })
    mat <- do.call(cbind, pcs)

    if (unitarize & !unitarized) {
      mat <- apply(mat, 2, .unitarize)
    } else if (unitarize & unitarized) {
      message("Already unitarized")
    }

    dt <- melt(
      as.data.table(mat)[, n := .I],
      id.vars = "n",
      variable.name = "name",
      value.name = "pc"
    )

    gr <- GRanges(rownames(mat))

    new_object(
      S7_object(),
      name = name,
      gr = gr,
      dt = dt,
      res = res,
      unitarized = unitarized,
      seqinfo = selectMethod('seqinfo', "GRanges")(gr),
      colnames = colnames(mat),
      mat = mat
    )
  }
)
S4_register(SingleCellCompartmentCall)

# Store the BiocGenerics %in% to differentiate from base::`%in%`
`%gin%` <- function(a, b) {
  BiocGenerics::match(a, b, nomatch = 0L) > 0L
}
