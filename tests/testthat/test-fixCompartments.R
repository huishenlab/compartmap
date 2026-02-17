min.conf <- 0.5

test_that(".inverter", {
  inverter <- compartmap:::.inverter
  row <- data.frame("conf.est" = c(1, 2, 3))
  expect_length(inverter(row, min.conf), nrow(row))
  expect_false(
    all(inverter(row, min.conf))
  )

  row <- 0.1 * row
  expect_true(
    all(inverter(row, min.conf))
  )
})

gr <- GRanges(c("chr1:1-10", "chr1:10-20", "chr1:21-30", "chr1:31-40", "chr1:41-50"))
conf.est <- c(0.1, 0.4, 0.5, 1, 2)
mcols(gr) <- list('score' = c(1, 2, 3, 4, 5))
mcols(gr)$conf.est <- conf.est

gr.no_conf_est <- gr
mcols(gr.no_conf_est) <- list(non_conf_est = conf.est)

gr1 <- gr
expected_mcols.gr1 <- DataFrame(
  score = c(1, 2, 3, 4, 5),
  conf.est = c(0.1, 0.4, 0.5, 1, 2),
  flip.compartment = as.logical(c(1, 1, 0, 0, 0)),
  flip.score = c(-1, -2, 3, 4, 5),
  flip.conf.est = c(0.9, 0.6, 0.5, 1, 2)
)

# score * 2
gr2 <- gr
mcols(gr2)$score <- mcols(gr)$score * 2
expected_mcols.gr2 <- DataFrame(
  score = c(2, 4, 6, 8, 10),
  conf.est = c(0.1, 0.4, 0.5, 1, 2),
  flip.compartment = as.logical(c(1, 1, 0, 0, 0)),
  flip.score = c(-2, -4, 6, 8, 10),
  flip.conf.est = c(0.9, 0.6, 0.5, 1, 2)
)

# conf.est * 2
gr3 <- gr
mcols(gr3)$conf.est <- mcols(gr)$conf.est * 2
expected_mcols.gr3 <- DataFrame(
  score = c(1, 2, 3, 4, 5),
  conf.est = c(0.2, 0.8, 1.0, 2.0, 4.0),
  flip.compartment = as.logical(c(1, 0, 0, 0, 0)),
  flip.score = c(-1, 2, 3, 4, 5),
  flip.conf.est = c(0.8, 0.8, 1.0, 2.0, 4.0)
)

gr1.expected <- gr
gr2.expected <- gr
gr3.expected <- gr

mcols(gr1.expected) <- expected_mcols.gr1
mcols(gr2.expected) <- expected_mcols.gr2
mcols(gr3.expected) <- expected_mcols.gr3

gr.input <- list(gr1, gr2, gr3)
gr.expected <- list(gr1.expected, gr2.expected, gr3.expected)

test_that("flipper", {
  expect_error(
    compartmap:::flipper(gr.no_conf_est, min.conf = 0.8),
    "Bootstrapping was not performed. Cannot fix compartments.",
    fixed = TRUE
  )

  lapply(1:length(gr.expected), function(i) {
    expect_equal(
      compartmap:::flipper(gr.input[[i]], min.conf),
      gr.expected[[i]]
    )
  })
})

test_that("fixCompartments", {
  lapply(1:length(gr.expected), function(i) {
    expect_equal(
      fixCompartments(gr.input[[i]], min.conf),
      gr.expected[[i]]
    )
  })

  re <- RaggedExperiment(A = gr1, B = gr2, C = gr3)
  fix_result <- fixCompartments(re)

  expect_equal(ncol(fix_result), ncol(re))
  expect_equal(colnames(fix_result), c("A", "B", "C"))
})
