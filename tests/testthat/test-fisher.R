test_that(".squeezit", {
  lapply(c(1:10), function(i) {
    expect_lt(compartmap:::.squeezeit(i), i)
    expect_equal(compartmap:::.squeezeit(i) / 0.999999, i)
  })
})

set.seed(42)
test_that("(i)fisherZ", {
  input <- matrix(sample(c(1, -1), 10000, replace = TRUE), ncol = 100)

  expect_equal(fisherZ(input), atanh(input * 0.999999))
  expect_equal(ifisherZ(input), tanh(input))

  expect_equal(fisherZ(ifisherZ(input)), input, tolerance = 0.00001)
  expect_equal(ifisherZ(fisherZ(input)), input, tolerance = 0.00001)
  expect_equal(fisherZ(ifisherZ(input)), input, tolerance = 0.00001)
})
