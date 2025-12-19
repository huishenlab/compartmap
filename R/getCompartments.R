#' Run compartment inference
#' @importFrom BiocParallel bplapply
#' @importFrom futile.logger flog.info flog.debug
getCompartments <- function(
  obj,
  res,
  chr,
  group,
  targets,
  bootstrap,
  num.bootstraps,
  genome,
  assay,
  boot.parallel,
  bpparams
) {
  if (is.null(chr)) {
    flog.info("Assuming we want to process all chromosomes.")
    chr <- getChrs(obj)
  }

  if (is.null(colnames(obj))) stop("colnames needs to be sample names.")
  columns <- colnames(obj)
  names(columns) <- columns

  prior.means <- getGlobalMeans(obj = obj, targets = targets, assay = assay)

  if (bootstrap) {
    flog.info("Pre-computing the bootstrap global means.")
    bmeans <- precomputeBootstrapMeans(
      obj = obj,
      BPPARAM = bpparams[[1]],
      targets = targets,
      num.bootstraps = num.bootstraps,
      assay = assay
    )
  }

  boot_msg <- ""
  if (bootstrap & boot.parallel) {
    boot_msg <- sprintf("Bootstrapping in parallel with %d cores", bpnworkers(bpparams[[2]]))
  } else if (bootstrap & !boot.parallel & group) {
    boot_msg <- "Not bootstrapping in parallel could take a long time..."
  }

  if (group) {
    flog.info("Computing group level compartments")
    flog.info(boot_msg)
    compartments.list <- bplapply(
      chr,
      function(c) {
        .getCompartments(
          obj,
          obj,
          assay = assay,
          BPPARAM = bpparams[[2]],
          res = res,
          chr = c,
          group = group,
          targets = targets,
          genome = genome,
          prior.means = prior.means,
          bootstrap = bootstrap,
          num.bootstraps = num.bootstraps,
          bootstrap.means = bmeans
        )
      },
      BPPARAM = bpparams[[1]]
    )

    compartments <- sort(unlist(as(compartments.list, "GRangesList")))
    return(compartments)
  }

  flog.info("Computing single-cell level compartments")
  flog.info(boot_msg)
  compartments <- bplapply(
    columns,
    function(s) {
      obj.sub <- obj[, s]

      compartments.list <- lapply(chr, function(c) {
      .getCompartments(
        obj.sub,
        obj,
        assay = assay,
        BPPARAM = bpparams[[2]],
        res = res,
        chr = c,
        group = group,
        targets = targets,
        genome = genome,
        prior.means = prior.means,
        bootstrap = bootstrap,
        num.bootstraps = num.bootstraps,
        bootstrap.means = bmeans
      )
    })
      sort(unlist(as(compartments.list, "GRangesList")))
    },
    BPPARAM = bpparams[[1]]
  )

  compartments <- as(compartments, "CompressedGRangesList")
  RaggedExperiment(compartments, colData = colData(obj))
}

# this is the main analysis function for computing compartments
.getCompartments <- function(
  obj,
  original.obj,
  assay,
  BPPARAM,
  res = 1e6,
  chr = NULL,
  group = FALSE,
  targets = NULL,
  genome = c("hg19", "hg38", "mm9", "mm10"),
  prior.means = NULL,
  bootstrap = TRUE,
  num.bootstraps = 1000,
  bootstrap.means = NULL
) {
  genome <- match.arg(genome)

  # update
  flog.debug("Computing compartments for %s", chr)
  obj <- keepSeqlevels(obj, chr, pruning.mode = "coarse")
  original.obj <- keepSeqlevels(original.obj, chr, pruning.mode = "coarse")

  # take care of the global means
  if (!is.null(prior.means)) {
    # this assumes that we've alread computed the global means
    pmeans <- as(prior.means, "GRanges")
    pmeans <- keepSeqlevels(pmeans, chr, pruning.mode = "coarse")
    # go back to a matrix
    prior.means <- as(pmeans, "matrix")
    colnames(prior.means) <- "globalMean"
  }

  obj.bins <- shrinkBins(
    obj,
    original.obj,
    prior.means = prior.means,
    chr = chr,
    res = res,
    targets = targets,
    assay = assay,
    genome = genome,
    jse = TRUE
  )

  obj.cor <- getCorMatrix(obj.bins, squeeze = !group)

  if (any(is.na(obj.cor$binmat.cor))) {
    obj.cor$gr$pc <- matrix(rep(NA, nrow(obj.cor$binmat.cor)))
    obj.svd <- obj.cor$gr
  } else {
    # compute SVD of correlation matrix
    obj.svd <- getABSignal(obj.cor, assay = assay)
  }

  if (!bootstrap) {
    return(obj.svd)
  }

  # bootstrap the estimates
  # always compute confidence intervals too
  # take care of the global means
  # this assumes that we've alread computed the global means
  bmeans <- as(bootstrap.means, "GRanges")
  bmeans <- keepSeqlevels(bmeans, chr, pruning.mode = "coarse")
  # go back to a matrix
  bmeans <- as(bmeans, "matrix")
  colnames(bmeans) <- rep("globalMean", ncol(bmeans))

  bootstrapCompartments(
    obj,
    original.obj,
    BPPARAM = BPPARAM,
    bootstrap.samples = num.bootstraps,
    chr = chr,
    group = group,
    assay = assay,
    targets = targets,
    res = res,
    genome = genome,
    q = 0.95,
    svd = obj.svd,
    bootstrap.means = bmeans
  )
}
