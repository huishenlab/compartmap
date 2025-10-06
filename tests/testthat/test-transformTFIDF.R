# transformTFIDF {{{
test_that("transformTFIDF", {
  expect_error(
    transformTFIDF(iris),
    "Input needs to be a matrix",
    fixed = TRUE
  )
  expect_error(
    transformTFIDF(0:10, count.min = 1, count.max = 0),
    "'count.min' must be less than 'count.max'"
  )
})
# }}}

# .constrain {{{
test_that(".constrain", {
  expect_equal(
    compartmap:::.constrain(0:10, lower = 0, upper = 1),
    c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  )
  expect_equal(
    compartmap:::.constrain(0:10, lower = 2, upper = 5),
    c(0, 0, 2, 3, 4, 5, 5, 5, 5, 5, 5)
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
