test_that("plotAB checks", {
  ranges <- c("chr1:1-10", "chr1:10-20")
  gr <- GRanges(ranges)
  cols <- list(
    'a' = c(1, 2),
    'b' = c(1, 2)
  )
  mcols(gr) <- cols
  expect_error(plotAB(matrix()), "'grAB' is not a GRanges object", fixed = TRUE)

  expect_error(
    plotAB(GRanges()),
    "score is not among names(mcols(x))",
    fixed = TRUE
  )
})


test_that(".unitarize", {
  vec <- c(1, 2, 3, 4, 5)
  vec.withNA <- c(1, 2, 3, NA, 5)

  expected <- {
    centered <- vec - median(vec, na.rm = TRUE)
    scaling <- sqrt(sum(centered^2))
    centered / scaling
  }
  expect_equal(compartmap:::.unitarize(vec), expected)

  expected.nomedian <- {
    scaling <- sqrt(sum(vec^2))
    vec / scaling
  }
  expect_equal(compartmap:::.unitarize(vec, medianCenter = FALSE), expected.nomedian)

  expected.withNA <- {
    centered <- vec.withNA - median(vec.withNA, na.rm = TRUE)
    scaling <- sqrt(sum(centered[-4]^2))
    centered / scaling
  }
  expect_message(
    compartmap:::.unitarize(vec.withNA, medianCenter = FALSE),
    "[.unitarize] 1 missing values were ignored.",
    fixed = TRUE
  )
  expect_equal(compartmap:::.unitarize(vec.withNA), expected.withNA)

  expected.withNA.nomedian <- {
    scaling <- sqrt(sum(vec.withNA[-4]^2, na.rm = TRUE))
    outvec <- vec.withNA
    outvec[c(1, 2, 3, 5)] <- vec.withNA[c(1, 2, 3, 5)] / scaling
    outvec
  }
  expect_equal(compartmap:::.unitarize(vec.withNA, medianCenter = FALSE), expected.withNA.nomedian)
})
