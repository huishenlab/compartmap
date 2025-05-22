test_that("getCorMatrix", {
  gr <- GRanges(c("chr1:1-10", "chr1:10-20", "chr1:21-30", "chr1:31-40"))
  count.mat <- matrix(rnbinom(16, 1, 0.5), ncol = 4)
  gmeans <- rowMeans(count.mat)
  binmat <- list(x = count.mat, gmeans = gmeans, gr = gr)

  expected.cormat <- t(binmat$x) |> cor()
  expected.cormat.squeezed <- compartmap:::fisherZ(expected.cormat)

  expected.result <- list(gr.cor = gr, binmat.cor = expected.cormat)
  expected.result.squeezed <- list(gr.cor = gr, binmat.cor = expected.cormat.squeezed)

  expect_message(getCorMatrix(binmat), "Calculating correlations")
  expect_message(getCorMatrix(binmat), "Done")
  expect_no_warning(expect_message(getCorMatrix(binmat)))

  expect_equal(getCorMatrix(binmat), expected.result)
  expect_equal(getCorMatrix(binmat, squeeze = TRUE), expected.result.squeezed)
})
