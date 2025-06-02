test_that("getSVD", {
})

library("BiocSingular")
test_that(".getUV", {
  mat <- matrix(runif(1000, 1, 100), nrow = 100)
  svd <- runSVD(mat, k = 1, BSPARAM = IrlbaParam())
  expect_equal(
    abs(compartmap:::.getUV(mat, k = 1, "left")),
    abs(svd$u)
  )
  expect_equal(
    abs(compartmap:::.getUV(mat, k = 1, "right")),
    abs(svd$v)
  )
})

test_that(".centerMatrix", {
  m <- matrix(rnorm(100), ncol = 10)

  expected <- t(apply(m, 1, function(i) {
    i - mean(i)
  }))

  expect_equal(
    compartmap:::.centerMatrix(m),
    expected,
    check.attributes = FALSE
  )
})
