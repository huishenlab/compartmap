
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

test_that(".getSummary", {
  gr1 <- GRanges(c("chr1:1-10", "chr1:11-20", "chr1:21-30", "chr1:31-40", "chr1:41-50", "chr1:51-60"))
  mcols(gr1) <- data.frame(score = runif(6, -5, 5))

  gr1$compartments <- ifelse(mcols(gr1)$score > 0, "open", "closed")
  gr1$open <- ifelse(mcols(gr1)$score > 0, 1, 0)
  gr1$closed <- ifelse(mcols(gr1)$score < 0, 1, 0)
  expected <- as.matrix(cbind(gr1$open, gr1$closed))

  expect_equal(compartmap:::.getSummary(gr1, gr1), expected)
})

test_that("summarizeBootstraps", {
})
