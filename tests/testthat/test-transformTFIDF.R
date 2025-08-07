# transformTFIDF {{{
test_that("transformTFIDF", {
  expect_error(
    transformTFIDF(iris),
    "Input needs to be a matrix",
    fixed = TRUE
  )
})
# }}}

# .tfidf {{{
test_that(".tfidf", {
  i <- c(1, 3:8)
  j <- c(2, 9, 6:10)
  x <- 1:7
  mat <- sparseMatrix(i, j, x = x) ##  8 x 10 "dgCMatrix"
  expect_equal(
    compartmap:::.tfidf(mat, 1:8),
    mat * 1:8
  )
})
# }}}
