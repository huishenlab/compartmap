test_that("getSVD", {})

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
