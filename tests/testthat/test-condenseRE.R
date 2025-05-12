grl <- GRangesList(
  GRanges(c("chr1:1-5", "chr1:4-6", "chr1:10-15"), score = 1:3, count = 1:3),
  GRanges(c("chr1:1-5", "chr2:1-3"), score = 4:5, count = 6:7)
)
names(grl) <- c("A", "B")
re <- RaggedExperiment(grl)
re.condense <- condenseRE(re)

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
