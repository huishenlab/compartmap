sample_names <- c("A", "B")
assay_names <- c("score", "count")
grl <- GRangesList(
  GRanges(c("chr1:1-5", "chr1:4-6", "chr1:10-15"), score = 1:3, count = 1:3),
  GRanges(c("chr1:1-5", "chr2:1-3"), score = 4:5, count = 6:7)
)
names(grl) <- sample_names
re <- RaggedExperiment(grl)
re.condenseRE <- condenseRE(re)
re.condenseSE <- condenseSE(re)

grl.single <- GRangesList(
  GRanges(c("chr1:1-5", "chr1:4-6", "chr1:10-15"), score = 1:3, count = 1:3)
)
names(grl.single) <- "A"
re.single <- RaggedExperiment(grl.single)
re.single.condenseSE <- condenseSE(re.single)

mat.score <- matrix(c(1, 2, 3, NA, 4, NA, NA, 5), nrow = 4)
rownames(mat.score) <- c("chr1:1-5", "chr1:4-6", "chr1:10-15", "chr2:1-3")
colnames(mat.score) <- c("A", "B")

mat.count <- matrix(c(1, 2, 3, NA, 6, NA, NA, 7), nrow = 4)
rownames(mat.count) <- c("chr1:1-5", "chr1:4-6", "chr1:10-15", "chr2:1-3")
colnames(mat.count) <- c("A", "B")

test_that("condenseRE", {
  obj <- SummarizedExperiment()
  expect_error(
    condenseRE(obj),
    "Input needs to be a RaggedExperiment"
  )

  res_list <- condenseRE(re)
  expect_equal(
    assayNames(re),
    names(res_list)
  )

  expect_equal(
    unname(rowRanges(res_list$score)),
    unique(rowRanges(re))
  )
  expect_equal(
    rowRanges(res_list$score),
    rowRanges(res_list$count)
  )

  expect_equal(
    assays(res_list$score)$score,
    mat.score
  )
  expect_equal(
    assays(res_list$count)$count,
    mat.count
  )
})

# condenseSE {{{
test_that("condenseSE", {
  obj <- RaggedExperiment()
  length(obj)
  expect_error(condenseSE(obj), "No assays found to condense.")

  expect_error(
    condenseSE(re, sample.name = "NA"),
    "sample.name <- sample.name %||% colnames(obj)"
  )

  expect_equal(length(re.condenseSE), 2)
  expect_equal(names(re.condenseSE), sample_names)

  lapply(1:length(re.condenseSE), function(i) {
    expected.gr <- unique(rowRanges(re))
    names(expected.gr) <- as.character(expected.gr)
    mcols(expected.gr)["score"] <- mat.score[, i]
    mcols(expected.gr)["count"] <- mat.count[, i]
    expect_equal(re.condenseSE[[i]], expected.gr)
  })

  expected.gr <- unique(rowRanges(re.single))
  names(expected.gr) <- as.character(expected.gr)
  mcols(expected.gr)["score"] <- mat.score[1:3, 1]
  mcols(expected.gr)["count"] <- mat.count[1:3, 1]
  expect_equal(re.single.condenseSE, expected.gr)
})
# }}}

