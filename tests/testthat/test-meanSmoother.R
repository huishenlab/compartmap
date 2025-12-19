test_that("meanSmoother", {
  mat <- matrix(1:5, nrow = 5, ncol = 5)
  expect_equal(
    meanSmoother(mat, k = 0),
    mat
  )

  # TODO get better error
  expect_error(
    meanSmoother(mat[1:2], k = 3),
    "length(mat) >= k is not TRUE",
    fixed = TRUE
  )

  # no NAs in output
  mat <- matrix(1:5, nrow = 5, ncol = 5)
  mat[5:6] <- NA
  expect_true(
    all(!is.na(meanSmoother(mat))),
  )
})

test_that(".window.mean", {
  mat <- matrix(1:16, ncol = 4)
  weights <- rep(1, length(mat))
  n <- length(mat)

  # k = 0 returns the matrix
  expect_equal(
    sapply(1:n, function(i) {
      compartmap:::.window.mean(mat, weights, i, k = 0)
    }),
    1:n
  )

  # k = 1 middle window
  sapply(2:(n - 1), function(i) {
    expect_equal(
      compartmap:::.window.mean(mat, weights, i, k = 1),
      mean(mat[(i - 1 - 1):(i + 1)])
    )
  })

  mat <- c(1:20)
  weights <- rep(1, length(mat))
  for (j in 3:10) {
    expect_equal(
      compartmap:::.window.mean(mat, weights, pos = j, k = 0),
      mat[j]
    )
    for (k in 2:j - 1) {
      expect_equal(
        compartmap:::.window.mean(mat, weights, pos = j, k = k),
        mean(mat[(j - k - 1):(j + k)])
      )
    }
  }
})

test_that(".window.mean edge cases", {
  mat <- c(1:20)
  weights <- rep(1, length(mat))

  # Test first position
  expect_equal(
    compartmap:::.window.mean(mat, weights, pos = 1, k = 0),
    mat[1]
  )

  # Test last position
  expect_equal(
    compartmap:::.window.mean(mat, weights, pos = 20, k = 0),
    mat[20]
  )

  # Test with NA values - should ignore NA
  mat_with_na <- c(1, 2, NA, 4, 5)
  expect_equal(
    compartmap:::.window.mean(mat_with_na, weights, pos = 3, k = 1),
    mean(c(1, 2, NA, 4), na.rm = TRUE)
  )

  # Test with different weights
  weights <- c(1, 2, 3, 4, 5)
  expect_equal(
    compartmap:::.window.mean(mat, weights, pos = 3, k = 1),
    weighted.mean(mat[1:4], weights[1:4]) # Weighted mean for position 3
  )
})

test_that(".meanSmoother.internal", {
  .meanSmoother.internal <- compartmap:::.meanSmoother.internal

  mat <- matrix(1:25, nrow = 5)
  weights.same <- rep(1, length(mat))
  weights.diff <- sample(1:2, size = length(mat), replace = TRUE)

  # k becomes 0 so should return 1
  window1 <- weighted.mean(mat[0:1], weights.same[0:1])

  # k remains 1
  window2 <- lapply(2:24, function(i) {
    stride <- (i - 1 - 1):(i + 1)
    weighted.mean(mat[stride], weights.same[stride])
  }) |>
    unlist()

  # k becomes 0 so should return 25
  window3 <- mat[25]

  expected_result <- c(window1, window2, window3)

  expect_equal(
    .meanSmoother.internal(mat, weights.same, k = 1),
    expected_result
  )

  window1 <- mat[1]
  window2 <- lapply(2:24, function(i) {
    stride <- (i - 1 - 1):(i + 1)
    weighted.mean(mat[stride], weights.diff[stride])
  }) |>
    unlist()
  window3 <- mat[25]

  expected_result <- c(window1, window2, window3)
  expect_equal(
    .meanSmoother.internal(mat, weights.diff, k = 1),
    expected_result
  )
})
