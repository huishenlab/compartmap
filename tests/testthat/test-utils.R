library(SummarizedExperiment)
library(RaggedExperiment)
library(GenomicRanges)

se <- SummarizedExperiment()
re <- RaggedExperiment()
df <- data.frame()
chrs <- c("chr2", "chr2", "chr1", "chr3")
gr <- GRanges(
  Rle(chrs, c(1, 3, 2, 4)),
  IRanges(1:10, width = 10:1)
)

# extractOpenClosed {{{
gr.pc <- gr
pc.list <- list('pc' = c(1, -1, 1, -1, -2, 2, -5, 5, -10, 10))
pc.state.cutoff.zero <- c(
  "open",
  "closed",
  "open",
  "closed",
  "closed",
  "open",
  "closed",
  "open",
  "closed",
  "open"
)
pc.state.cutoff.two <- c(
  "closed",
  "closed",
  "closed",
  "closed",
  "closed",
  "closed",
  "closed",
  "open",
  "closed",
  "open"
)
pc.state.array <- lapply(
  pc.state.cutoff.zero,
  function(i) ifelse(i == "open", "closed", "open")
) |>
  unlist()
mcols(gr.pc) <- pc.list

test_that("extractOpenClosed", {
  expect_error(extractOpenClosed(se))
  expect_error(extractOpenClosed(re))
  expect_error(extractOpenClosed(df))
  expect_error(extractOpenClosed(gr))
  expect_no_error(extractOpenClosed(gr.pc))
  expect_equal(
    extractOpenClosed(gr.pc, assay = "rna"),
    pc.state.cutoff.zero
  )
  expect_equal(
    extractOpenClosed(gr.pc, assay = "rna", cutoff = 2),
    pc.state.cutoff.two
  )
  expect_equal(
    extractOpenClosed(gr.pc, assay = "atac"),
    pc.state.cutoff.zero
  )
  expect_equal(
    extractOpenClosed(gr.pc, assay = "array"),
    pc.state.array
  )
})
# }}}

# verifySE{{{
test_that("verifySE", {
  expect_no_error(compartmap:::verifySE(se))
  expect_error(compartmap:::verifySE(re), err.verifySE)
  expect_error(compartmap:::verifySE(gr), err.verifySE)
  expect_error(compartmap:::verifySE(df), err.verifySE)
})
# }}}

# verifyCoords {{{
test_that("verifyCoords", {
  se.noranges <- SummarizedExperiment(rowRanges = GRanges())
  se.ranges <- SummarizedExperiment(rowRanges = gr)
  expect_error(
    compartmap:::verifyCoords(se.noranges),
    err.verifyCoords
  )
  expect_no_error(compartmap:::verifyCoords(se.ranges))
})
# }}}

# verifyAssayNames{{{
se.rna <- SummarizedExperiment(
  assays = SimpleList(counts = matrix()),
)
se.atac <- SummarizedExperiment(
  assays = SimpleList(counts = matrix()),
)
se.array <- SummarizedExperiment(
  assays = SimpleList(Beta = matrix()),
)
se.bisulfite <- SummarizedExperiment(
  assays = SimpleList(Beta = matrix()),
)

test_that("verifyAssayNames", {
  expect_error(
    compartmap:::verifyAssayNames(se, "s"),
    paste(shQuote('s'), "is unsupported")
  )
  expect_error(
    compartmap:::verifyAssayNames(se.rna, "array"),
    err.verifyAssayNames.beta
  )
  expect_error(
    compartmap:::verifyAssayNames(se.array, "rna"),
    err.verifyAssayNames.rna_counts
  )
  expect_no_error(
    compartmap:::verifyAssayNames(se.rna, "rna")
  )
  expect_no_error(
    compartmap:::verifyAssayNames(se.atac, "atac")
  )
  expect_no_error(
    compartmap:::verifyAssayNames(se.array, "array")
  )
  expect_no_error(
    compartmap:::verifyAssayNames(se.bisulfite, "bisulfite")
  )
})
# }}}

# getAssayNames{{{

nrows <- 200
ncols <- 6
counts1 <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
counts2 <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
counts3 <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
colData <- DataFrame(
  Treatment = rep(c("a", "b"), 3),
  row.names = LETTERS[1:6]
)
se <- SummarizedExperiment(
  assays = SimpleList(one = counts1, two = counts2, three = counts3),
  colData = colData
)

test_that("getAssayNames", {
  expect_equal(
    getAssayNames(se),
    c("one", "two", "three")
  )
})
# }}}

# flogit {{{
test_that("flogit", {
  input <- c(1, 0, 0.5, -0.5, -1, 1.01)
  expect_error(flogit(input), "Values of p outside (0,1) detected.", fixed = T)
})
# }}}

# fexpit {{{
test_that("fexpit", {
  input <- rnorm(10)
  output <- fexpit(input)
  expect_true(all(output > 0))
  expect_true(all(output < 1))
})
# }}}

# getChrs {{{
test_that("getChrs", {
  expect_equal(
    getChrs(gr),
    unique(chrs)
  )
})
# }}}

# removeEmptyBoots {{{
test_that("removeEmptyBoots", {
  gr1 <- GRanges("chr2", IRanges(3, 6), strand = "+", score = 5L, GC = 0.45)
  gr2 <- GRanges(
    c("chr1", "chr1"),
    IRanges(c(7, 13), width = 3),
    strand = c("+", "-"),
    score = 3:4,
    GC = c(0.3, 0.5)
  )
  gr3 <- GRanges(
    c("chr1", "chr2"),
    IRanges(c(1, 4), c(3, 9)),
    strand = c("-", "-"),
    score = c(6L, 2L),
    GC = c(0.4, 0.1)
  )
  grlist <- list(gr1, gr2, gr3)
  grlist.na <- list(gr1, gr2, NA, gr3, NA)

  expect_equal(
    removeEmptyBoots(grlist),
    grlist
  )
  expect_equal(
    removeEmptyBoots(grlist.na),
    grlist
  )
})
# }}}

# getGenome {{{
test_that("getGenome", {
  bundled_genomes <- c("hg19", "hg38", "mm9", "mm10")
  lapply(bundled_genomes, function(g) {
    expect_equal(
      getGenome(g),
      get(paste0(g, ".gr"))
    )
    expect_equal(
      getGenome(g, type = "openseas"),
      get(paste0("openSeas.", g))
    )
  })
  expect_error(getGenome("NA50"))
})
# }}}

# getSeqLengths{{{
forseqlengths.gr <- GRanges(
  Rle(c("chr1", "chr2", "chr3"), c(1, 3, 6)),
  IRanges(1:10, width = 10:1, names = head(letters, 10)),
  Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  score = 1:10,
  GC = seq(1, 0, length = 10)
)
seqlengths.list <- c(chr1 = 100, chr2 = 200, chr3 = 300)
seqinfo(forseqlengths.gr) <- GenomeInfoDb::Seqinfo(
  paste0("chr", 1:3),
  seqlengths.list,
  NA,
  "mock1"
)
# to get seqlengths(gr) set to
# chr1 chr2 chr3
#  100  200  300

seqlengths.list <- c(chr1 = 100, chr2 = 200, chr3 = 300)
seqinfo(forseqlengths.gr) <- GenomeInfoDb::Seqinfo(
  paste0("chr", 1:3),
  seqlengths.list,
  NA,
  "mock1"
)
test_that("getSeqLengths", {
  lapply(seq(1:3), function(i) {
    expect_equal(
      getSeqLengths(forseqlengths.gr, paste0("chr", i)),
      seqlengths.list[i]
    )
  })
  # missing chromosome
  expect_error(getSeqLengths(forseqlengths.gr, "chr0"))
})
# }}}

# cleanAssay {{{

checkCleanAssay <- function(se.array, se.bisulfite, cleanFun, dimFun) {
  expect_equal(
    dimFun(assays(cleanFun(se.array, assay = "array"))$Beta),
    0
  )
  expect_equal(
    dimFun(assays(cleanFun(se.bisulfite, assay = "bisulfite"))$counts),
    0
  )
  expect_equal(
    dimFun(assays(cleanFun(se.array, assay = "array", na.max = 0.9))$Beta),
    10
  )
  expect_equal(
    dimFun(
      assays(cleanFun(se.bisulfite, assay = "bisulfite", na.max = 0.9))$counts
    ),
    10
  )
}

test_that("cleanAssayRows", {
  expect_type(
    compartmap:::cleanAssay(by = "row"),
    "closure"
  )
  expect_type(
    compartmap:::cleanAssay(by = "col"),
    "closure"
  )
})

test_that("cleanAssayRows", {
  cleanFun <- compartmap:::cleanAssayRows
  dimFun <- nrow
  m <- matrix(1:10, 10, 10)
  m.allNA <- apply(m, 1, function(i) ifelse(i <= 2, i, NA))
  se.array <- SummarizedExperiment(assays = SimpleList(Beta = m.allNA))
  se.bisulfite <- SummarizedExperiment(assays = SimpleList(counts = m.allNA))
  checkCleanAssay(se.array, se.bisulfite, cleanFun, dimFun)
})

test_that("cleanAssayCols", {
  cleanFun <- compartmap:::cleanAssayCols
  dimFun <- ncol
  m <- matrix(1:10, 10, 10)
  m.allNA <- apply(m, 2, function(i) ifelse(i <= 2, i, NA))
  se.array <- SummarizedExperiment(assays = SimpleList(Beta = m.allNA))
  se.bisulfite <- SummarizedExperiment(assays = SimpleList(counts = m.allNA))
  checkCleanAssay(se.array, se.bisulfite, cleanFun, dimFun)
})
### }}}

# filterOpenSea {{{
test_that("filterOpenSea", {
  expect_error(
    filterOpenSea(df, genome = "hg38"),
    "'obj' needs to be a GRanges or SummarizedExperiment",
    fixed = TRUE
  )
  expect_error(
    filterOpenSea(gr, other = se, genome = "hg38"),
    "The 'other' input needs to be a GRanges of open sea regions",
    fixed = TRUE
  )

  nrows <- 10
  ncols <- 6
  counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
  colData <- DataFrame(
    Treatment = rep(c("ChIP", "Input"), 3),
    row.names = LETTERS[1:6]
  )
  se <- SummarizedExperiment(
    assays = SimpleList(counts = counts),
    colData = colData,
    rowRanges = gr
  )

  rownames.cg <- paste0(rep("cg", 5), 1:5)
  rownames.non_cg <- 6:10
  rownames(se) <- c(rownames.cg, rownames.non_cg)
  expected.se <- se[1:5, ]

  expect_equal(
    filterOpenSea(se, genome = "hg19"),
    expected.se
  )
})
# }}}
