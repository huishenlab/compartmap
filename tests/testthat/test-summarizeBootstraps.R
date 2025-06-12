
test_that(".getCI", {
  gr <- GRanges(c(
    "chr1:1-10",
    "chr1:10-20",
    "chr1:21-30",
    "chr1:31-40",
    "chr1:41-50",
    "chr1:51-60"
  ))

  est.ab <- data.frame(
    compartments = rep(c("open", "closed"), 3),
    boot.open = c(100, 50, 5, 2, 1, 0),
    boot.closed = c(0, 50, 95, 98, 99, 100)
  )
  mcols(gr) <- est.ab

  gr.expected <- gr
  successes <- c(100, 50, 5, 98, 1, 100)
  failures <- c(0, 50, 95, 2, 99, 0)
  expectedCI <- agrestiCoullCI(successes, failures, 0.95)
  mcols(gr.expected) <- cbind(
    mcols(gr),
    agrestiCoullCI(successes, failures, 0.95)
  )

  expect_equal(gr.expected, compartmap:::.getCI(gr, 0.95))
})


test_that("summarizeBootstraps", {
})
