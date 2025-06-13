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
  # create the non-bootstrapped calls as all open
  est.ab <- GRanges(c(
    "chr1:1-10",
    "chr1:11-20",
    "chr1:21-30",
    "chr1:31-40",
    "chr1:41-50",
    "chr1:51-60"
  ))
  mcols(est.ab) <- data.frame(pc = rep(1, 6))
  est.ab$compartments <- ifelse(mcols(est.ab)$pc > 0, "open", "closed")
  gr1 <- est.ab
  gr2 <- est.ab
  mcols(gr2[3, ]) <- NA

  removeEmptyBoots(list(gr1, gr2))

  expect_message(summarizeBootstraps(list(gr1, gr1), gr1), "Summarizing bootstraps")

  tester <- function(one, two, expected) {
    mcols(gr1) <- data.frame(pc = one)
    gr1$compartments <- ifelse(mcols(gr1)$pc > 0, "open", "closed")

    mcols(gr2) <- data.frame(pc = two)
    gr2$compartments <- ifelse(mcols(gr2)$pc > 0, "open", "closed")

    boot.list <- list(gr1, gr2)
    expect_equal(
      mcols(summarizeBootstraps(boot.list, est.ab))[, 6:8],
      expected
    )
  }

  # set up the successes and failures based on the all open est.ab
  test_set <- list(
    # equal successes and failures
    list(one = rep(-1, 6), two = rep(1, 6), ci = agrestiCoullCI(1, 1, 0.95)),
    # only successes
    list(one = rep(1, 6), two = rep(1, 6), ci = agrestiCoullCI(2, 0, 0.95)),
    # only failures
    list(one = rep(0, 6), two = rep(0, 6), ci = agrestiCoullCI(0, 2, 0.95))
  )

  lapply(test_set, function(set) {
    ci <- set$ci
    expected <- rbind(ci, ci, ci, ci, ci, ci) |> DataFrame()
    tester(set$one, set$two, expected)
  })

  # alternating success/failures on one
  mcols(gr1) <- data.frame(pc = rep(c(-1, 1), 3))
  gr1$compartments <- ifelse(mcols(gr1)$pc > 0, "open", "closed")

  mcols(gr2) <- data.frame(pc = rep(0, 6))
  gr2$compartments <- ifelse(mcols(gr2)$pc > 0, "open", "closed")

  boot.list <- list(gr1, gr2)

  expected <- rbind(
    agrestiCoullCI(0, 2, 0.95),
    agrestiCoullCI(1, 1, 0.95),
    agrestiCoullCI(0, 2, 0.95),
    agrestiCoullCI(1, 1, 0.95),
    agrestiCoullCI(0, 2, 0.95),
    agrestiCoullCI(1, 1, 0.95)
  ) |>
    DataFrame()

  expect_equal(
    mcols(summarizeBootstraps(boot.list, est.ab))[, 6:8],
    expected
  )
})
