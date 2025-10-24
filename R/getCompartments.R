getCompartments <- function(
  obj,
  assay,
  res,
  parallel,
  chr,
  targets,
  cores,
  bootstrap,
  num.bootstraps,
  boot.parallel,
  boot.cores,
  genome,
  group
) {
  if (is.null(chr)) {
    message("Assuming we want to process all chromosomes.")
    chr <- getChrs(obj)
  }

  if (is.null(colnames(obj))) stop("colnames needs to be sample names.")
  columns <- colnames(obj)
  names(columns) <- columns

  prior.means <- getGlobalMeans(obj = obj, targets = targets, assay = assay)

  if (bootstrap) {
    message("Pre-computing the bootstrap global means.")
    bmeans <- precomputeBootstrapMeans(
      obj = obj,
      targets = targets,
      num.bootstraps = num.bootstraps,
      assay = assay,
      parallel = parallel,
      num.cores = cores
    )
  }

  if (group) {
    compartments.list <- mclapply(chr, function(c) {
      .getCompartments(
        obj,
        obj,
        assay = assay,
        res = res,
        chr = c,
        targets = targets,
        genome = genome,
        bootstrap = bootstrap,
        num.bootstraps = num.bootstraps,
        prior.means = prior.means,
        parallel = boot.parallel,
        cores = boot.cores,
        group = group,
        bootstrap.means = bmeans
      )
    }, mc.cores = ifelse(parallel, cores, 1))

    compartments <- sort(unlist(as(compartments.list, "GRangesList")))
    return(compartments)
  }

  compartments <- mclapply(columns, function(s) {
    obj.sub <- obj[, s]

    message("Working on ", s)
    compartments.list <- lapply(chr, function(c) {
      .getCompartments(
        obj.sub,
        obj,
        assay = assay,
        res = res,
        chr = c,
        targets = targets,
        genome = genome,
        bootstrap = bootstrap,
        prior.means = prior.means,
        num.bootstraps = num.bootstraps,
        parallel = boot.parallel,
        cores = boot.cores,
        group = group,
        bootstrap.means = bmeans
      )
    })
    sort(unlist(as(compartments.list, "GRangesList")))
  }, mc.cores = ifelse(parallel, cores, 1), mc.preschedule = F)

  compartments <- as(compartments, "CompressedGRangesList")
  RaggedExperiment(compartments, colData = colData(obj))
}

# this is the main analysis function for computing compartments
.getCompartments <- function(
  obj,
  original.obj,
  assay,
  res = 1e6,
  chr = NULL,
  targets = NULL,
  genome = c("hg19", "hg38", "mm9", "mm10"),
  prior.means = NULL,
  bootstrap = TRUE,
  num.bootstraps = 1000,
  parallel = FALSE,
  cores = 2,
  group = FALSE,
  bootstrap.means = NULL
) {
  genome <- match.arg(genome)

  if (parallel) options(mc.cores = cores)

  # update
  message("Computing compartments for ", chr)
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
    bootstrap.samples = num.bootstraps,
    chr = chr,
    assay = assay,
    parallel = parallel,
    cores = cores,
    targets = targets,
    res = res,
    genome = genome,
    q = 0.95,
    svd = obj.svd,
    group = group,
    bootstrap.means = bmeans
  )
}
