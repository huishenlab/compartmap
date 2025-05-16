test_that("flipper", {
  gr <- GRanges(c("chr1:1-10", "chr1:10-20"))
  cols <- list('b' = c(1, 2))
  mcols(gr) <- cols
  expect_error(
    compartmap:::flipper(gr, min.conf = 0.8),
    "Bootstrapping was not performed. Cannot fix compartments.",
    fixed = TRUE
  )

  gr.conf <- gr
  cols.conf <- list('conf.est' = c(1, 2))
  mcols(gr.conf) <- cols.conf
})

test_that(".inverter", {
  inverter <- compartmap:::.inverter
  row <- data.frame("conf.est" = c(1, 2, 3))
  min.conf <- 0.5
  expect_length(inverter(row, min.conf), nrow(row))
  expect_false(
    all(inverter(row, min.conf))
  )

  row <- 0.1 * row
  expect_true(
    all(inverter(row, min.conf))
  )
})
